#!/usr/bin/env python3
"""
Standalone tangle detection from scaffold GAF files.
Parses scaffolds for [N...N:*] markers, finds boundary nodes by walking
outward until nodes meet length/coverage thresholds, clusters gaps across
haplotypes, and outputs detected tangles.
"""

import sys
import os
import re
import json
import logging
import argparse
import networkx as nx
from dataclasses import dataclass, field
from collections import defaultdict, deque
from typing import Optional

# Add project root to path so we can import src modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.node_id_mapper import NodeIdMapper
from src.input_parsing import parse_gfa, read_coverage_file, coverage_from_graph, verify_coverage, get_nonoriented_graph, get_nonoriented_dual_graph
from src.graph_transformation import create_dual_graph

# Constants
MIN_BOUNDARY_LENGTH = 40_000       # Minimum node length to qualify as boundary
COV_LOW_FACTOR      = 0.5          # Boundary cov >= median * COV_LOW_FACTOR
COV_HIGH_FACTOR     = 1.5          # Boundary cov <= median * COV_HIGH_FACTOR
GAP_PATTERN         = re.compile(r'^\[N\d+N:.*\]$')
MIN_SCAFFOLD_LENGTH = 5_000_000    # 5 Mbp: minimum scaffold length for multichromosomal check
MAX_WALK_DISTANCE   = 5_000_000      # Max cumulative node length when searching for boundaries


@dataclass
class ScaffoldToken:
    """A single token from scaffold path: either an oriented node or a gap marker."""
    raw: str
    is_gap: bool = False
    orientation: str = ''
    node_name: str = ''


@dataclass
class GapInfo:
    """One gap found in a single scaffold line."""
    scaffold_name: str
    gap_marker: str
    gap_index: int = 0
    left_boundary: Optional[str] = None
    right_boundary: Optional[str] = None
    left_orientation: str = ''
    right_orientation: str = ''
    inner_node_names: list = field(default_factory=list)


@dataclass
class DetectedTangle:
    """A tangle detected from one or more gaps across haplotypes."""
    tangle_id: int
    gaps: list = field(default_factory=list)
    boundary_pairs: list = field(default_factory=list)
    inner_nodes: set = field(default_factory=set)
    is_multihaplotype: bool = False
    is_multichromosomal: bool = False
    has_no_graph_path: bool = False
    has_shared_boundary: bool = False
    boundaries_do_not_separate: bool = False
    boundaries_not_synchronized: bool = False
    haplotypes: set = field(default_factory=set)
    all_scaffolds: set = field(default_factory=set)
    notes: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core utilities
# ---------------------------------------------------------------------------

def tokenize_scaffold_path(path_str):
    """Parse a scaffold path string into tokens.

    Supports two formats:
      GAF:  <node1>node2[N100N:scaffold]<node3   (orientation prefix, no separators)
      TSV:  node1-,node2+,[N100N:scaffold],node3- (orientation suffix +/-, comma-separated)
    """
    # Detect format: comma-separated means TSV style
    if ',' in path_str:
        return _tokenize_csv_path(path_str)
    return _tokenize_gaf_path(path_str)


def _tokenize_gaf_path(path_str):
    """Parse GAF-style path: <node>node[gap]."""
    tokens = []
    i = 0
    n = len(path_str)
    while i < n:
        if path_str[i] == '[':
            j = path_str.index(']', i)
            raw = path_str[i:j+1]
            tokens.append(ScaffoldToken(raw=raw, is_gap=True))
            i = j + 1
        elif path_str[i] in '<>':
            orient = path_str[i]
            j = i + 1
            while j < n and path_str[j] not in '<>[':
                j += 1
            node_name = path_str[i+1:j]
            tokens.append(ScaffoldToken(
                raw=path_str[i:j], is_gap=False,
                orientation=orient, node_name=node_name))
            i = j
        else:
            i += 1
    return tokens


def _tokenize_csv_path(path_str):
    """Parse CSV-style path: node1+,node2-,[N100N:type],node3+."""
    tokens = []
    for part in path_str.split(','):
        part = part.strip()
        if not part:
            continue
        if part.startswith('['):
            tokens.append(ScaffoldToken(raw=part, is_gap=True))
        elif part.endswith('+') or part.endswith('-'):
            orient = '>' if part.endswith('+') else '<'
            node_name = part[:-1]
            tokens.append(ScaffoldToken(
                raw=part, is_gap=False,
                orientation=orient, node_name=node_name))
        else:
            # Unoriented node — treat as forward
            tokens.append(ScaffoldToken(
                raw=part, is_gap=False,
                orientation='>', node_name=part))
    return tokens


def is_real_gap(token):
    """Check if token is a [N...N:*] marker."""
    return token.is_gap and bool(GAP_PATTERN.match(token.raw))


def node_qualifies_as_boundary(node_name, graph, cov, median_cov, node_mapper):
    """Check if a node meets length + coverage thresholds."""
    if not node_mapper.has_name(node_name):
        return False
    node_id = node_mapper.get_id_for_name(node_name)
    if node_id not in graph.nodes:
        return False
    length = graph.nodes[node_id].get('length', 0)
    coverage = cov.get(node_id, 0)
    return (length >= MIN_BOUNDARY_LENGTH and
            median_cov * COV_LOW_FACTOR <= coverage <= median_cov * COV_HIGH_FACTOR)


def is_in_graph(node_name, node_mapper, graph):
    """Check if node_name exists in graph."""
    if not node_mapper.has_name(node_name):
        return False
    node_id = node_mapper.get_id_for_name(node_name)
    return node_id in graph.nodes or (-node_id) in graph.nodes


def get_haplotype(scaffold_name):
    """Extract haplotype identifier from scaffold name."""
    parts = scaffold_name.split('_from_')
    return parts[0] if parts else scaffold_name


def oriented_node_id(node_name, orientation, node_mapper):
    """Get the signed node id from name + orientation."""
    if not node_mapper.has_name(node_name):
        return None
    nid = node_mapper.get_id_for_name(node_name)
    return -nid if orientation == '<' else nid


# ---------------------------------------------------------------------------
# Graph path analysis
# ---------------------------------------------------------------------------

def find_graph_path(left_name, left_orient, right_name, right_orient,
                    graph, node_mapper, max_depth=200):
    """
    BFS from left boundary to right boundary in directed graph.
    Returns (path_exists, inner_node_ids) where inner_node_ids is a set of
    absolute node IDs on the shortest path (excluding the boundaries themselves).
    """
    start_id = oriented_node_id(left_name, left_orient, node_mapper)
    end_id = oriented_node_id(right_name, right_orient, node_mapper)
    if start_id is None or end_id is None:
        return False, set()
    if start_id not in graph.nodes or end_id not in graph.nodes:
        return False, set()

    visited = {}
    queue = deque()
    for succ in graph.successors(start_id):
        queue.append((succ, 1))
        visited[succ] = start_id

    while queue:
        node, depth = queue.popleft()
        if node == end_id:
            # Trace back path
            path_nodes = set()
            cur = node
            while cur != start_id:
                path_nodes.add(abs(cur))
                cur = visited[cur]
            path_nodes.discard(abs(start_id))
            path_nodes.discard(abs(end_id))
            return True, path_nodes
        if depth >= max_depth:
            continue
        for succ in graph.successors(node):
            if succ not in visited:
                visited[succ] = node
                queue.append((succ, depth + 1))

    return False, set()


def find_inner_nodes_from_graph(boundary_pairs, seed_inner_names, graph, nonoriented_dual, node_mapper):
    """
    Find all graph nodes enclosed between boundary nodes using the dual graph.
    In the dual graph, original nodes are edges. We remove edges corresponding
    to boundary nodes, then find the connected component containing a seed
    inner node's edge. This correctly handles node-based GFA graphs where the
    tangle interior may be disconnected in the original graph but connected
    through shared junctions in the dual graph.
    Returns set of node names.
    """
    boundary_ids = set()
    for bp in boundary_pairs:
        for key in ('start', 'end'):
            name = bp[key]
            if node_mapper.has_name(name):
                boundary_ids.add(node_mapper.get_id_for_name(name))

    if not boundary_ids:
        return seed_inner_names

    # Remove dual edges corresponding to boundary nodes
    subgraph = nonoriented_dual.copy()
    edges_to_remove = []
    for u, v, key, data in subgraph.edges(keys=True, data=True):
        orig = abs(data.get('original_node', 0))
        if orig in boundary_ids:
            edges_to_remove.append((u, v, key))
    for u, v, key in edges_to_remove:
        subgraph.remove_edge(u, v, key=key)

    # Find a seed junction vertex from known inner nodes
    seed_vertex = None
    for name in seed_inner_names:
        if node_mapper.has_name(name):
            nid = node_mapper.get_id_for_name(name)
            # Find a dual edge for this inner node
            for u, v, data in nonoriented_dual.edges(data=True):
                if abs(data.get('original_node', 0)) == nid:
                    # Use endpoint that's still connected in subgraph
                    if u in subgraph and subgraph.degree(u) > 0:
                        seed_vertex = u
                        break
                    if v in subgraph and subgraph.degree(v) > 0:
                        seed_vertex = v
                        break
            if seed_vertex is not None:
                break

    # If no seed from walk, try BFS path nodes (using directed graph)
    if seed_vertex is None:
        for bp in boundary_pairs:
            _, path_ids = find_graph_path(
                bp['start'], bp['start_orientation'],
                bp['end'], bp['end_orientation'],
                graph, node_mapper)
            for pid in path_ids:
                if pid in boundary_ids:
                    continue
                for u, v, data in nonoriented_dual.edges(data=True):
                    if abs(data.get('original_node', 0)) == pid:
                        if u in subgraph and subgraph.degree(u) > 0:
                            seed_vertex = u
                            break
                        if v in subgraph and subgraph.degree(v) > 0:
                            seed_vertex = v
                            break
                if seed_vertex is not None:
                    break
            if seed_vertex is not None:
                break

    # If no seed, try boundary neighbor junctions (smallest component)
    if seed_vertex is None:
        # Collect junction vertices adjacent to boundary edges
        boundary_junctions = set()
        for u, v, key in edges_to_remove:
            boundary_junctions.add(u)
            boundary_junctions.add(v)
        seen = set()
        best_component = None
        for j in boundary_junctions:
            if j in seen or j not in subgraph or subgraph.degree(j) == 0:
                continue
            comp = nx.node_connected_component(subgraph, j)
            seen.update(comp)
            if best_component is None or len(comp) < len(best_component):
                best_component = comp
        if best_component is not None:
            result = set()
            for u, v, data in subgraph.edges(data=True):
                if u in best_component and v in best_component:
                    nid = abs(data['original_node'])
                    name = node_mapper.node_id_to_name_safe(nid)
                    if name:
                        result.add(name.lstrip('<>'))
            return result
        return seed_inner_names

    if seed_vertex not in subgraph or subgraph.degree(seed_vertex) == 0:
        return seed_inner_names

    component = nx.node_connected_component(subgraph, seed_vertex)
    result = set()
    for u, v, data in subgraph.edges(data=True):
        if u in component and v in component:
            nid = abs(data['original_node'])
            if nid not in boundary_ids:
                name = node_mapper.node_id_to_name_safe(nid)
                if name:
                    result.add(name.lstrip('<>'))
    return result


def validate_tangle_boundaries(boundary_pairs, inner_node_names, graph, nonoriented_dual, node_mapper):
    """
    Validate that boundary nodes properly separate the tangle from the rest
    of the graph. Uses the dual graph where boundary nodes are edges:
    removing boundary edges should disconnect the inner component from the
    outer component. Also checks that each boundary node separates
    incoming/outgoing traffic in the directed graph.
    Returns (is_valid, notes).
    """
    notes = []

    boundary_node_ids = set()
    boundary_orientations = {}
    for bp in boundary_pairs:
        for key, orient_key in [('start', 'start_orientation'), ('end', 'end_orientation')]:
            name = bp[key]
            if node_mapper.has_name(name):
                nid = node_mapper.get_id_for_name(name)
                boundary_node_ids.add(nid)
                orient = bp[orient_key]
                oriented_id = -nid if orient == '<' else nid
                if nid not in boundary_orientations:
                    boundary_orientations[nid] = set()
                boundary_orientations[nid].add(oriented_id)

    inner_node_ids = set()
    for name in inner_node_names:
        if node_mapper.has_name(name):
            inner_node_ids.add(node_mapper.get_id_for_name(name))

    if not boundary_node_ids or not inner_node_ids:
        return True, notes

    # Remove dual edges corresponding to boundary nodes
    subgraph = nonoriented_dual.copy()
    boundary_edges = []
    for u, v, key, data in subgraph.edges(keys=True, data=True):
        orig = abs(data.get('original_node', 0))
        if orig in boundary_node_ids:
            boundary_edges.append((u, v, key))
    for u, v, key in boundary_edges:
        subgraph.remove_edge(u, v, key=key)

    # Find a seed junction vertex for an inner node
    seed_vertex = None
    for nid in inner_node_ids:
        for u, v, data in nonoriented_dual.edges(data=True):
            if abs(data.get('original_node', 0)) == nid:
                if u in subgraph and subgraph.degree(u) > 0:
                    seed_vertex = u
                    break
                if v in subgraph and subgraph.degree(v) > 0:
                    seed_vertex = v
                    break
        if seed_vertex is not None:
            break

    if seed_vertex is None:
        notes.append("No inner nodes found in dual graph")
        return False, notes

    tangle_component = nx.node_connected_component(subgraph, seed_vertex)

    # Collect inner node IDs from tangle component edges
    tangle_node_ids = set()
    for u, v, data in subgraph.edges(data=True):
        if u in tangle_component and v in tangle_component:
            tangle_node_ids.add(abs(data['original_node']))

    is_valid = True

    # Check directional separation: each oriented boundary should have
    # predecessors and successors on DIFFERENT sides of the boundary.
    for bnid, oriented_ids in boundary_orientations.items():
        for oriented_id in oriented_ids:
            if oriented_id not in graph.nodes:
                continue
            has_pred_in = any(abs(p) in tangle_node_ids for p in graph.predecessors(oriented_id))
            has_succ_in = any(abs(s) in tangle_node_ids for s in graph.successors(oriented_id))
            if has_pred_in and has_succ_in:
                bname = node_mapper.node_id_to_name_safe(bnid)
                orient_char = '>' if oriented_id > 0 else '<'
                notes.append(
                    f"Boundary {orient_char}{bname} does not separate incoming/outgoing")
                is_valid = False

    return is_valid, notes


# ---------------------------------------------------------------------------
# Scaffold analysis
# ---------------------------------------------------------------------------

def walk_for_boundary(tokens, gap_idx, direction, graph, cov, median_cov, node_mapper):
    """
    Walk from gap_idx in given direction (+1=right, -1=left).
    Returns (boundary_node_name, orientation, inner_nodes) or (None, '', inner_nodes).
    """
    inner_nodes = []
    idx = gap_idx + direction
    while 0 <= idx < len(tokens):
        t = tokens[idx]
        if is_real_gap(t):
            return None, '', inner_nodes
        if not t.is_gap:
            if is_in_graph(t.node_name, node_mapper, graph):
                if node_qualifies_as_boundary(t.node_name, graph, cov, median_cov, node_mapper):
                    return t.node_name, t.orientation, inner_nodes
                else:
                    inner_nodes.append(t.node_name)
        idx += direction
    return None, '', inner_nodes


def detect_gaps_in_scaffold(scaffold_name, path_str, graph, cov, median_cov, node_mapper):
    """Find all gap markers in one scaffold and walk outward for boundaries."""
    tokens = tokenize_scaffold_path(path_str)
    gap_indices = [i for i, t in enumerate(tokens) if is_real_gap(t)]
    if not gap_indices:
        return []

    gaps = []
    for gap_num, gi in enumerate(gap_indices, start=1):
        gap = GapInfo(scaffold_name=scaffold_name, gap_marker=tokens[gi].raw, gap_index=gap_num)

        left_name, left_orient, left_inner = walk_for_boundary(
            tokens, gi, -1, graph, cov, median_cov, node_mapper)
        gap.left_boundary = left_name
        gap.left_orientation = left_orient

        right_name, right_orient, right_inner = walk_for_boundary(
            tokens, gi, +1, graph, cov, median_cov, node_mapper)
        gap.right_boundary = right_name
        gap.right_orientation = right_orient

        gap.inner_node_names = left_inner + right_inner
        gaps.append(gap)

    return gaps


def cluster_gaps_into_tangles(all_gaps, graph, nonoriented_graph, nonoriented_dual, node_mapper):
    """Cluster gaps into tangles using union-find on shared inner nodes."""
    n = len(all_gaps)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    node_to_gaps = defaultdict(list)
    for i, gap in enumerate(all_gaps):
        # Only use inner nodes for clustering if a graph path exists
        # between the boundaries; otherwise the walked nodes aren't
        # truly enclosed in a tangle.
        if gap.left_boundary and gap.right_boundary:
            path_exists, _ = find_graph_path(
                gap.left_boundary, gap.left_orientation,
                gap.right_boundary, gap.right_orientation,
                graph, node_mapper)
            if not path_exists:
                continue
        for name in gap.inner_node_names:
            node_to_gaps[name].append(i)

    for name, gap_indices in node_to_gaps.items():
        for j in range(1, len(gap_indices)):
            union(gap_indices[0], gap_indices[j])

    # Union consecutive gaps in the same scaffold when no proper boundary
    # separates them: either a boundary is missing, or consecutive gaps
    # share the same boundary node (one's right_boundary == next's left_boundary).
    scaffold_gap_indices = defaultdict(list)
    for i, gap in enumerate(all_gaps):
        scaffold_gap_indices[gap.scaffold_name].append(i)

    for scaffold_name, indices in scaffold_gap_indices.items():
        indices.sort(key=lambda idx: all_gaps[idx].gap_index)
        for k in range(len(indices) - 1):
            gi = indices[k]
            gj = indices[k + 1]
            gap_i = all_gaps[gi]
            gap_j = all_gaps[gj]
            # Merge only when both boundaries are missing (no boundary
            # exists between consecutive gaps in the scaffold).
            if gap_i.right_boundary is None and gap_j.left_boundary is None:
                union(gi, gj)

    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)

    tangles = []
    for tid, (root, indices) in enumerate(sorted(clusters.items())):
        cluster_gaps = [all_gaps[i] for i in indices]

        inner = set()
        for g in cluster_gaps:
            inner.update(g.inner_node_names)

        haplotypes = set()
        for g in cluster_gaps:
            haplotypes.add(get_haplotype(g.scaffold_name))

        notes = []

        is_multihaplotype = len(haplotypes) > 2
        if is_multihaplotype:
            notes.append(f"Multihaplotype ({len(haplotypes)} haplotypes: {', '.join(sorted(haplotypes))}): unresolvable")

        # Check for shared boundaries
        has_shared_boundary = False
        boundary_to_haplotypes = defaultdict(set)
        for g in cluster_gaps:
            hap = get_haplotype(g.scaffold_name)
            if g.left_boundary:
                boundary_to_haplotypes[g.left_boundary].add(hap)
            if g.right_boundary:
                boundary_to_haplotypes[g.right_boundary].add(hap)
        shared_nodes = [node for node, haps in boundary_to_haplotypes.items() if len(haps) > 1]
        if shared_nodes:
            has_shared_boundary = True
            for sn in shared_nodes:
                haps = sorted(boundary_to_haplotypes[sn])
                notes.append(f"Shared boundary node {sn} used by haplotypes: {', '.join(haps)}")

        # Build boundary pairs per scaffold by aggregating: use the leftmost
        # boundary of the first gap and rightmost boundary of the last gap.
        # For single-gap scaffolds this is equivalent to taking the gap's pair.
        # For multi-gap scaffolds (merged adjacent gaps) this yields the
        # correct outer boundaries.
        seen_pairs = set()
        boundary_pairs = []
        scaffold_gaps_map = defaultdict(list)
        for g in cluster_gaps:
            scaffold_gaps_map[g.scaffold_name].append(g)

        for scaffold_name, s_gaps in scaffold_gaps_map.items():
            left_boundary = left_orient = right_boundary = right_orient = None
            for g in s_gaps:
                if g.left_boundary:
                    left_boundary = g.left_boundary
                    left_orient = g.left_orientation
                    break
            for g in reversed(s_gaps):
                if g.right_boundary:
                    right_boundary = g.right_boundary
                    right_orient = g.right_orientation
                    break

            if left_boundary and right_boundary:
                pair_key = (left_boundary, right_boundary)
                if pair_key not in seen_pairs:
                    seen_pairs.add(pair_key)
                    boundary_pairs.append({
                        'start': left_boundary,
                        'start_orientation': left_orient,
                        'end': right_boundary,
                        'end_orientation': right_orient,
                        'scaffold': scaffold_name,
                    })
                    if len(s_gaps) > 1:
                        logging.info(
                            f"  Tangle cluster: aggregated boundaries from {len(s_gaps)} "
                            f"gaps in {scaffold_name}: "
                            f"{left_orient}{left_boundary} -> {right_orient}{right_boundary}")

        # Validate boundaries
        boundaries_do_not_separate = False
        if boundary_pairs:
            is_valid, validation_notes = validate_tangle_boundaries(
                boundary_pairs, inner, graph, nonoriented_dual, node_mapper)
            if not is_valid:
                boundaries_do_not_separate = True
            notes.extend(validation_notes)

        # Check graph path existence
        has_no_graph_path = False
        for bp in boundary_pairs:
            path_exists, _ = find_graph_path(
                bp['start'], bp['start_orientation'],
                bp['end'], bp['end_orientation'],
                graph, node_mapper)
            if not path_exists:
                has_no_graph_path = True
                notes.append(
                    f"No graph path from {bp['start_orientation']}{bp['start']} to "
                    f"{bp['end_orientation']}{bp['end']} "
                    f"(scaffold {bp.get('scaffold', '?')})")

        if not boundary_pairs:
            for g in cluster_gaps:
                if g.left_boundary is None and g.right_boundary is None:
                    notes.append(f"Both boundaries missing (scaffold {g.scaffold_name})")

        # Compute inner nodes from graph
        if boundary_pairs:
            inner = find_inner_nodes_from_graph(boundary_pairs, inner, graph, nonoriented_dual, node_mapper)

        tangles.append(DetectedTangle(
            tangle_id=tid,
            gaps=cluster_gaps,
            boundary_pairs=boundary_pairs,
            inner_nodes=inner,
            is_multihaplotype=is_multihaplotype,
            has_no_graph_path=has_no_graph_path,
            has_shared_boundary=has_shared_boundary,
            boundaries_do_not_separate=boundaries_do_not_separate,
            haplotypes=haplotypes,
            notes=notes,
        ))

    return tangles


def compute_median_coverage(graph, cov):
    """Length-weighted median coverage over all positive-id nodes."""
    nodes = []
    total_length = 0
    for node_id in graph.nodes():
        if node_id < 0:
            continue
        length = graph.nodes[node_id].get('length', 0)
        coverage = cov.get(node_id, 0)
        nodes.append((length, coverage))
        total_length += length

    nodes.sort(key=lambda x: x[1])
    cum = 0
    for length, coverage in nodes:
        cum += length
        if cum * 2 >= total_length:
            return coverage
    return 0.0


def resolve_shared_boundaries(tangle, all_scaffold_tokens, graph, cov, median_cov, node_mapper):
    """
    When a boundary node is shared across haplotypes, walk further outward
    in each scaffold to find per-haplotype distinct boundaries.
    """
    if not tangle.has_shared_boundary:
        return

    boundary_to_haplotypes = defaultdict(set)
    boundary_to_pairs = defaultdict(list)
    for i, bp in enumerate(tangle.boundary_pairs):
        hap = get_haplotype(bp['scaffold'])
        boundary_to_haplotypes[bp['start']].add(hap)
        boundary_to_haplotypes[bp['end']].add(hap)
        boundary_to_pairs[bp['start']].append((i, 'start'))
        boundary_to_pairs[bp['end']].append((i, 'end'))

    shared_nodes = {node for node, haps in boundary_to_haplotypes.items() if len(haps) > 1}
    if not shared_nodes:
        return

    for shared_node in shared_nodes:
        logging.info(f"  Tangle {tangle.tangle_id}: resolving shared boundary {shared_node}")
        tangle.inner_nodes.add(shared_node)

        for pair_idx, side in boundary_to_pairs[shared_node]:
            bp = tangle.boundary_pairs[pair_idx]
            scaffold_name = bp['scaffold']
            if scaffold_name not in all_scaffold_tokens:
                continue

            tokens = all_scaffold_tokens[scaffold_name]
            shared_token_indices = [
                i for i, t in enumerate(tokens)
                if not t.is_gap and t.node_name == shared_node
            ]
            if not shared_token_indices:
                continue

            direction = -1 if side == 'start' else +1
            shared_idx = shared_token_indices[0] if side == 'start' else shared_token_indices[-1]

            new_boundary = None
            new_orient = ''
            extra_inner = []
            idx = shared_idx + direction
            while 0 <= idx < len(tokens):
                t = tokens[idx]
                if is_real_gap(t):
                    break
                if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                    if node_qualifies_as_boundary(t.node_name, graph, cov, median_cov, node_mapper):
                        new_boundary = t.node_name
                        new_orient = t.orientation
                        break
                    else:
                        extra_inner.append(t.node_name)
                idx += direction

            if new_boundary:
                logging.info(
                    f"    {scaffold_name}: {side} boundary "
                    f"{bp.get(side + '_orientation', '')}{shared_node} -> "
                    f"{new_orient}{new_boundary}")
                bp[side] = new_boundary
                bp[side + '_orientation'] = new_orient
                tangle.inner_nodes.update(extra_inner)
            else:
                logging.warning(
                    f"    {scaffold_name}: could not find replacement for "
                    f"shared {side} boundary {shared_node}")

    # Re-check shared boundaries after resolution
    boundary_to_haplotypes = defaultdict(set)
    for bp in tangle.boundary_pairs:
        hap = get_haplotype(bp['scaffold'])
        boundary_to_haplotypes[bp['start']].add(hap)
        boundary_to_haplotypes[bp['end']].add(hap)
    still_shared = [n for n, h in boundary_to_haplotypes.items() if len(h) > 1]
    tangle.has_shared_boundary = bool(still_shared)
    tangle.notes = [n for n in tangle.notes if not n.startswith("Shared boundary")]
    for sn in still_shared:
        haps = sorted(boundary_to_haplotypes[sn])
        tangle.notes.append(f"Shared boundary node {sn} used by haplotypes: {', '.join(haps)}")


# ---------------------------------------------------------------------------
# Tangle classification
# ---------------------------------------------------------------------------

def extend_failing_boundaries(tangle, all_scaffold_tokens, graph, nonoriented_graph, nonoriented_dual,
                              cov, median_cov, node_mapper, max_rounds=5):
    """
    When boundary nodes fail directional separation, walk further outward
    in the scaffold to find replacement boundaries that properly separate.
    The old boundary becomes an inner node of the (now larger) tangle.
    Stops when cumulative absorbed length exceeds MAX_WALK_DISTANCE.
    """
    # Track cumulative absorbed length per (pair_idx, side)
    absorbed_length = defaultdict(int)

    for round_num in range(max_rounds):
        if not tangle.boundaries_do_not_separate or not tangle.boundary_pairs:
            return

        # Build tangle component by removing boundary edges from dual graph
        boundary_ids = set()
        for bp in tangle.boundary_pairs:
            for key in ('start', 'end'):
                name = bp[key]
                if node_mapper.has_name(name):
                    boundary_ids.add(node_mapper.get_id_for_name(name))

        inner_ids = set()
        for name in tangle.inner_nodes:
            if node_mapper.has_name(name):
                inner_ids.add(node_mapper.get_id_for_name(name))

        if not boundary_ids or not inner_ids:
            return

        # Remove boundary edges from dual graph
        sub_dual = nonoriented_dual.copy()
        b_edges = []
        for u, v, key, data in sub_dual.edges(keys=True, data=True):
            if abs(data.get('original_node', 0)) in boundary_ids:
                b_edges.append((u, v, key))
        for u, v, key in b_edges:
            sub_dual.remove_edge(u, v, key=key)

        # Find seed junction vertex for an inner node
        seed_vertex = None
        for nid in inner_ids:
            for u, v, data in nonoriented_dual.edges(data=True):
                if abs(data.get('original_node', 0)) == nid:
                    if u in sub_dual and sub_dual.degree(u) > 0:
                        seed_vertex = u
                        break
                    if v in sub_dual and sub_dual.degree(v) > 0:
                        seed_vertex = v
                        break
            if seed_vertex is not None:
                break
        if seed_vertex is None:
            return

        dual_component = nx.node_connected_component(sub_dual, seed_vertex)
        # Collect inner node IDs from dual component edges
        tangle_component = set()
        for u, v, data in sub_dual.edges(data=True):
            if u in dual_component and v in dual_component:
                tangle_component.add(abs(data['original_node']))

        # Find which (pair_idx, side) fail directional check
        failing = set()
        for pair_idx, bp in enumerate(tangle.boundary_pairs):
            for side, orient_key in [('start', 'start_orientation'),
                                     ('end', 'end_orientation')]:
                name = bp[side]
                orient = bp[orient_key]
                if not node_mapper.has_name(name):
                    continue
                nid = node_mapper.get_id_for_name(name)
                oriented_id = -nid if orient == '<' else nid
                if oriented_id not in graph.nodes:
                    continue
                has_pred_in = any(abs(p) in tangle_component
                                  for p in graph.predecessors(oriented_id))
                has_succ_in = any(abs(s) in tangle_component
                                  for s in graph.successors(oriented_id))
                if has_pred_in and has_succ_in:
                    failing.add((pair_idx, side))

        if not failing:
            return

        extended_any = False
        for pair_idx, side in failing:
            bp = tangle.boundary_pairs[pair_idx]
            old_boundary = bp[side]

            # Check if absorbing the old boundary would exceed distance limit.
            # Always allow the first extension attempt (absorbed_length==0)
            # even if the boundary node itself is very large.
            if node_mapper.has_name(old_boundary):
                old_nid = node_mapper.get_id_for_name(old_boundary)
                old_len = graph.nodes.get(old_nid, {}).get('length', 0)
            else:
                old_len = 0
            if (absorbed_length[(pair_idx, side)] > 0 and
                    absorbed_length[(pair_idx, side)] + old_len > MAX_WALK_DISTANCE):
                continue

            scaffold_name = bp['scaffold']

            if scaffold_name not in all_scaffold_tokens:
                continue

            tokens = all_scaffold_tokens[scaffold_name]
            boundary_token_indices = [
                i for i, t in enumerate(tokens)
                if not t.is_gap and t.node_name == old_boundary
            ]
            if not boundary_token_indices:
                continue

            direction = -1 if side == 'start' else +1
            start_idx = (boundary_token_indices[0] if side == 'start'
                         else boundary_token_indices[-1])

            new_boundary = None
            new_orient = ''
            extra_inner = []
            walk_length = 0
            idx = start_idx + direction
            while 0 <= idx < len(tokens):
                t = tokens[idx]
                if is_real_gap(t):
                    break
                if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                    nid = node_mapper.get_id_for_name(t.node_name)
                    node_len = graph.nodes.get(nid, {}).get('length', 0)
                    if node_qualifies_as_boundary(t.node_name, graph, cov,
                                                  median_cov, node_mapper):
                        new_boundary = t.node_name
                        new_orient = t.orientation
                        break
                    else:
                        walk_length += node_len
                        if absorbed_length[(pair_idx, side)] + walk_length > MAX_WALK_DISTANCE:
                            break
                        extra_inner.append(t.node_name)
                idx += direction

            if new_boundary:
                absorbed_length[(pair_idx, side)] += old_len + walk_length
                logging.info(
                    f"  Tangle {tangle.tangle_id}: extending {side} boundary "
                    f"{bp[side + '_orientation']}{old_boundary} -> "
                    f"{new_orient}{new_boundary} (scaffold {scaffold_name})")
                bp[side] = new_boundary
                bp[side + '_orientation'] = new_orient
                tangle.inner_nodes.add(old_boundary)
                tangle.inner_nodes.update(extra_inner)
                extended_any = True

        if not extended_any:
            return

        # Recompute inner nodes from graph with new boundaries
        tangle.inner_nodes = find_inner_nodes_from_graph(
            tangle.boundary_pairs, tangle.inner_nodes,
            graph, nonoriented_dual, node_mapper)

        # Re-validate
        tangle.notes = [n for n in tangle.notes
                        if not n.startswith("Boundaries do not isolate")
                        and "does not separate incoming/outgoing" not in n]
        is_valid, validation_notes = validate_tangle_boundaries(
            tangle.boundary_pairs, tangle.inner_nodes,
            graph, nonoriented_dual, node_mapper)
        tangle.boundaries_do_not_separate = not is_valid
        tangle.notes.extend(validation_notes)


def classify_tangle(t):
    """
    Classify a tangle into a category and compute flags.
    Returns (category, flags, is_invalid) where:
      category: 'valid_1hap', 'valid_2hap', 'no_path', 'multiscaffold', 'other_invalid'
      flags: list of flag strings
      is_invalid: bool
    """
    flags = []
    if t.is_multichromosomal:
        flags.append("MULTICHROMOSOMAL")
    if t.is_multihaplotype:
        flags.append("MULTIHAPLOTYPE/UNRESOLVABLE")
    if t.boundaries_not_synchronized:
        flags.append("NOT SYNCHRONIZED")
    if t.has_no_graph_path:
        flags.append("NO GRAPH PATH")
    if t.has_shared_boundary:
        flags.append("SHARED BOUNDARY")

    has_boundaries = bool(t.boundary_pairs)
    is_other_invalid = (not t.has_no_graph_path and
                        not t.is_multichromosomal and
                        not t.is_multihaplotype and
                        (t.boundaries_do_not_separate or
                         t.boundaries_not_synchronized or
                         t.has_shared_boundary or
                         not has_boundaries))
    if is_other_invalid:
        flags.append("OTHER INVALID")

    is_invalid = (t.has_no_graph_path or t.is_multichromosomal or
                  t.is_multihaplotype or is_other_invalid)

    if t.has_no_graph_path:
        category = 'no_path'
    elif t.is_multichromosomal or t.is_multihaplotype:
        category = 'multiscaffold'
    elif is_other_invalid:
        category = 'other_invalid'
    elif len(t.haplotypes) <= 1:
        category = 'valid_1hap'
    else:
        category = 'valid_2hap'

    return category, flags, is_invalid


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def detect_tangles_from_scaffolds(scaffold_file, graph, cov, median_cov, node_mapper):
    """Full pipeline: parse scaffolds, find gaps, walk boundaries, cluster."""

    nonoriented_graph = get_nonoriented_graph(graph)

    # Create dual graph (vertices=junctions, edges=original nodes) for
    # connectivity checks.  In the dual representation boundary nodes become
    # edges; removing them properly disconnects the tangle interior even when
    # the interior is not connected in the original node-based GFA.
    dual_graph = create_dual_graph(graph, node_mapper)
    nonoriented_dual = get_nonoriented_dual_graph(dual_graph, node_mapper)

    # Step 1: Parse and tokenize all scaffolds
    all_scaffold_tokens = {}
    all_scaffold_paths = {}
    with open(scaffold_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            scaffold_name = parts[0]
            path_str = parts[1]
            if scaffold_name == 'name':
                continue
            all_scaffold_tokens[scaffold_name] = tokenize_scaffold_path(path_str)
            all_scaffold_paths[scaffold_name] = path_str

    # Step 2: Detect gaps in each scaffold
    all_gaps = []
    for scaffold_name, path_str in all_scaffold_paths.items():
        gaps = detect_gaps_in_scaffold(
            scaffold_name, path_str, graph, cov, median_cov, node_mapper)
        all_gaps.extend(gaps)

    logging.info(f"Found {len(all_gaps)} gaps across all scaffolds")
    logging.info(f"Total scaffolds processed: {len(all_scaffold_paths)}")

    # Step 3: Cluster gaps into tangles
    tangles = cluster_gaps_into_tangles(all_gaps, graph, nonoriented_graph, nonoriented_dual, node_mapper)
    logging.info(f"Clustered into {len(tangles)} tangles")

    # Step 3b: Merge tangles that share a significant fraction of inner nodes.
    # This catches cases where gaps from unrelated scaffolds land in the same
    # graph region but weren't merged by union-find (e.g. because no graph
    # path existed between their boundaries at clustering time).
    # Tangles with no directed graph path between boundaries are excluded:
    # their "inner nodes" are just walked nodes, not a true tangle component.
    logging.info("Step 3b: Merging tangles with overlapping inner nodes...")
    merged_any = True
    while merged_any:
        merged_any = False
        for i in range(len(tangles)):
            if tangles[i] is None:
                continue
            ti = tangles[i]
            if not ti.inner_nodes or ti.has_no_graph_path:
                continue
            for j in range(i + 1, len(tangles)):
                if tangles[j] is None:
                    continue
                tj = tangles[j]
                if not tj.inner_nodes or tj.has_no_graph_path:
                    continue
                overlap = ti.inner_nodes & tj.inner_nodes
                larger = max(len(ti.inner_nodes), len(tj.inner_nodes))
                # Inner node sets should either coincide or be disjoint.
                # Merge only when they approximately coincide (>80% of
                # the larger set overlaps).
                if larger >= 10 and len(overlap) >= larger * 0.8:
                    # Merge tj into ti
                    logging.info(
                        f"  Merging tangles {ti.tangle_id} and {tj.tangle_id}: "
                        f"{len(overlap)} shared inner nodes "
                        f"(|ti|={len(ti.inner_nodes)}, |tj|={len(tj.inner_nodes)})")
                    all_cluster_gaps = ti.gaps + tj.gaps
                    # Build new boundary pairs per scaffold
                    scaffold_gaps_map = defaultdict(list)
                    for g in all_cluster_gaps:
                        scaffold_gaps_map[g.scaffold_name].append(g)
                    new_boundary_pairs = []
                    seen_bp = set()
                    for sn, s_gaps in scaffold_gaps_map.items():
                        left_b = left_o = right_b = right_o = None
                        for g in s_gaps:
                            if g.left_boundary:
                                left_b = g.left_boundary
                                left_o = g.left_orientation
                                break
                        for g in reversed(s_gaps):
                            if g.right_boundary:
                                right_b = g.right_boundary
                                right_o = g.right_orientation
                                break
                        if left_b and right_b:
                            pkey = (left_b, right_b)
                            if pkey not in seen_bp:
                                seen_bp.add(pkey)
                                new_boundary_pairs.append({
                                    'start': left_b, 'start_orientation': left_o,
                                    'end': right_b, 'end_orientation': right_o,
                                    'scaffold': sn,
                                })
                    new_inner = ti.inner_nodes | tj.inner_nodes
                    if new_boundary_pairs:
                        new_inner = find_inner_nodes_from_graph(
                            new_boundary_pairs, new_inner, graph,
                            nonoriented_dual, node_mapper)
                    new_haplotypes = ti.haplotypes | tj.haplotypes
                    # Re-check properties
                    new_has_no_path = False
                    new_notes = []
                    for bp in new_boundary_pairs:
                        path_exists, _ = find_graph_path(
                            bp['start'], bp['start_orientation'],
                            bp['end'], bp['end_orientation'],
                            graph, node_mapper)
                        if not path_exists:
                            new_has_no_path = True
                            new_notes.append(
                                f"No graph path from {bp['start_orientation']}{bp['start']} to "
                                f"{bp['end_orientation']}{bp['end']} "
                                f"(scaffold {bp.get('scaffold', '?')})")
                    bdns = False
                    has_shared = False
                    if new_boundary_pairs and new_inner:
                        is_valid, v_notes = validate_tangle_boundaries(
                            new_boundary_pairs, new_inner, graph,
                            nonoriented_dual, node_mapper)
                        bdns = not is_valid
                        new_notes.extend(v_notes)
                    boundary_to_haps = defaultdict(set)
                    for g in all_cluster_gaps:
                        hap = get_haplotype(g.scaffold_name)
                        if g.left_boundary:
                            boundary_to_haps[g.left_boundary].add(hap)
                        if g.right_boundary:
                            boundary_to_haps[g.right_boundary].add(hap)
                    if any(len(h) > 1 for h in boundary_to_haps.values()):
                        has_shared = True
                    tangles[i] = DetectedTangle(
                        tangle_id=ti.tangle_id,
                        gaps=all_cluster_gaps,
                        boundary_pairs=new_boundary_pairs,
                        inner_nodes=new_inner,
                        is_multihaplotype=len(new_haplotypes) > 2,
                        has_no_graph_path=new_has_no_path,
                        has_shared_boundary=has_shared,
                        boundaries_do_not_separate=bdns,
                        haplotypes=new_haplotypes,
                        notes=new_notes,
                    )
                    tangles[j] = None
                    merged_any = True
                    break
            if merged_any:
                break
        if merged_any:
            tangles = [t for t in tangles if t is not None]
            for idx, t in enumerate(tangles):
                t.tangle_id = idx
    logging.info(f"After inner-node merge: {len(tangles)} tangles")

    # Step 4: Resolve shared boundaries by walking further in scaffolds
    logging.info("Step 4: Resolving shared boundaries...")
    for t in tangles:
        if t.has_shared_boundary:
            resolve_shared_boundaries(t, all_scaffold_tokens, graph, cov, median_cov, node_mapper)

    # Step 5: Flag multichromosomal tangles
    logging.info("Step 5: Checking for multichromosomal tangles...")

    scaffold_length = {}
    for scaffold_name, tokens in all_scaffold_tokens.items():
        if '_from_' not in scaffold_name:
            continue
        total = 0
        for t_tok in tokens:
            if not t_tok.is_gap and t_tok.node_name and node_mapper.has_name(t_tok.node_name):
                nid = node_mapper.get_id_for_name(t_tok.node_name)
                for sid in (nid, -nid):
                    if sid in graph.nodes:
                        total += graph.nodes[sid].get('length', 0)
                        break
        scaffold_length[scaffold_name] = total

    node_to_scaffolds = defaultdict(set)
    for scaffold_name, tokens in all_scaffold_tokens.items():
        if '_from_' not in scaffold_name:
            continue
        if scaffold_length.get(scaffold_name, 0) < MIN_SCAFFOLD_LENGTH:
            continue
        for t_tok in tokens:
            if not t_tok.is_gap and t_tok.node_name:
                node_to_scaffolds[t_tok.node_name].add(scaffold_name)

    for t in tangles:
        gap_scaffolds = set(g.scaffold_name for g in t.gaps)
        t.all_scaffolds = set(gap_scaffolds)
        t.passthrough_scaffolds = set()
        for bp in t.boundary_pairs:
            t.all_scaffolds.add(bp['scaffold'])

        long_gap_scaffolds = [s for s in gap_scaffolds
                              if scaffold_length.get(s, 0) >= MIN_SCAFFOLD_LENGTH]
        if not long_gap_scaffolds:
            continue

        # If the tangle itself has gaps from >2 scaffolds, it's
        # multichromosomal — multiple chromosomes converge in this region.
        if len(gap_scaffolds) > 2:
            t.is_multichromosomal = True
            t.notes.append(
                f"Multichromosomal: gaps from {len(gap_scaffolds)} scaffolds: "
                f"{', '.join(sorted(gap_scaffolds))}")
            logging.info(
                f"  Tangle {t.tangle_id}: multichromosomal "
                f"({len(gap_scaffolds)} gap scaffolds)")
            continue

        passthrough_scaffolds = set()
        for inner_name in t.inner_nodes:
            for scaff in node_to_scaffolds.get(inner_name, set()):
                if scaff not in gap_scaffolds:
                    passthrough_scaffolds.add(scaff)

        t.all_scaffolds.update(passthrough_scaffolds)
        t.passthrough_scaffolds = passthrough_scaffolds
        if passthrough_scaffolds:
            total_scaffolds = gap_scaffolds | passthrough_scaffolds
            if len(total_scaffolds) > 2:
                t.is_multichromosomal = True
                t.notes.append(
                    f"Multichromosomal: scaffolds passing through without a gap: "
                    f"{', '.join(sorted(passthrough_scaffolds))}")
                logging.info(
                    f"  Tangle {t.tangle_id}: multichromosomal "
                    f"(passthrough: {', '.join(sorted(passthrough_scaffolds))})")
            else:
                logging.info(
                    f"  Tangle {t.tangle_id}: passthrough from "
                    f"{', '.join(sorted(passthrough_scaffolds))} "
                    f"(only {len(total_scaffolds)} scaffolds, not multichromosomal)")

    # Step 5.5: Add boundaries from passthrough scaffolds
    # When boundaries don't separate and a passthrough scaffold exists, the
    # second haplotype traverses the tangle region without a gap.  Find
    # where it diverges from the gap scaffold's path and add boundary pairs.
    logging.info("Step 5.5: Adding boundaries from passthrough scaffolds...")
    for t in tangles:
        if not t.boundaries_do_not_separate or not t.boundary_pairs:
            continue
        if t.is_multichromosomal:
            continue
        gap_scaffolds = set(g.scaffold_name for g in t.gaps)
        pt_scaffolds = getattr(t, 'passthrough_scaffolds', set())
        if not pt_scaffolds:
            continue

        added_any = False
        for pt_scaffold in sorted(pt_scaffolds):
            if pt_scaffold not in all_scaffold_tokens:
                continue
            pt_tokens = all_scaffold_tokens[pt_scaffold]

            for bp in list(t.boundary_pairs):
                # Only consider boundary pairs from gap scaffolds
                if bp['scaffold'] not in gap_scaffolds:
                    continue
                scaffold_name = bp['scaffold']
                if scaffold_name not in all_scaffold_tokens:
                    continue
                gap_tokens = all_scaffold_tokens[scaffold_name]

                # Find boundary positions in gap scaffold
                start_idx = None
                end_idx = None
                for i, tok in enumerate(gap_tokens):
                    if not tok.is_gap and tok.node_name == bp['start']:
                        start_idx = i
                    if not tok.is_gap and tok.node_name == bp['end']:
                        end_idx = i

                if start_idx is None or end_idx is None:
                    continue

                # Find anchor nodes: first token outside each boundary
                # on the gap scaffold (away from the tangle interior).
                pre_start = None
                idx = start_idx - 1 if start_idx < end_idx else start_idx + 1
                step = -1 if start_idx < end_idx else +1
                while 0 <= idx < len(gap_tokens):
                    tok = gap_tokens[idx]
                    if tok.is_gap:
                        break
                    if tok.node_name:
                        pre_start = tok.node_name
                        break
                    idx += step

                post_end = None
                idx = end_idx + 1 if start_idx < end_idx else end_idx - 1
                step = +1 if start_idx < end_idx else -1
                while 0 <= idx < len(gap_tokens):
                    tok = gap_tokens[idx]
                    if tok.is_gap:
                        break
                    if tok.node_name:
                        post_end = tok.node_name
                        break
                    idx += step

                if not pre_start or not post_end:
                    continue

                # Locate anchors on passthrough scaffold
                pt_pre_idx = None
                pt_post_idx = None
                for i, tok in enumerate(pt_tokens):
                    if not tok.is_gap and tok.node_name == pre_start:
                        pt_pre_idx = i
                    if not tok.is_gap and tok.node_name == post_end:
                        pt_post_idx = i

                if pt_pre_idx is None or pt_post_idx is None:
                    continue

                # Scan between anchors for boundary-qualifying nodes
                lo = min(pt_pre_idx, pt_post_idx)
                hi = max(pt_pre_idx, pt_post_idx)
                candidates = []
                for i in range(lo + 1, hi):
                    tok = pt_tokens[i]
                    if tok.is_gap or not tok.node_name:
                        continue
                    if not is_in_graph(tok.node_name, node_mapper, graph):
                        continue
                    if node_qualifies_as_boundary(tok.node_name, graph, cov,
                                                  median_cov, node_mapper):
                        candidates.append((tok.node_name, tok.orientation))

                if len(candidates) < 2:
                    continue

                # Assign start/end based on direction
                if pt_pre_idx < pt_post_idx:
                    pt_start_name, pt_start_orient = candidates[0]
                    pt_end_name, pt_end_orient = candidates[-1]
                else:
                    pt_start_name, pt_start_orient = candidates[-1]
                    pt_end_name, pt_end_orient = candidates[0]

                # Skip if these nodes are already boundaries
                existing = set()
                for ebp in t.boundary_pairs:
                    existing.add(ebp['start'])
                    existing.add(ebp['end'])
                if pt_start_name in existing and pt_end_name in existing:
                    continue

                new_bp = {
                    'start': pt_start_name,
                    'start_orientation': pt_start_orient,
                    'end': pt_end_name,
                    'end_orientation': pt_end_orient,
                    'scaffold': pt_scaffold,
                }
                t.boundary_pairs.append(new_bp)
                t.haplotypes.add(get_haplotype(pt_scaffold))
                added_any = True
                logging.info(
                    f"  Tangle {t.tangle_id}: added passthrough boundary pair "
                    f"{pt_start_orient}{pt_start_name} -> "
                    f"{pt_end_orient}{pt_end_name} "
                    f"(from scaffold {pt_scaffold})")

        if added_any:
            # Use original gap inner nodes as seeds (not the old expanded set)
            # to ensure the seed lands in the correct tangle component.
            seed_inner = set()
            for g in t.gaps:
                seed_inner.update(g.inner_node_names)
            # Recompute inner nodes with expanded boundary set
            t.inner_nodes = find_inner_nodes_from_graph(
                t.boundary_pairs, seed_inner if seed_inner else t.inner_nodes,
                graph, nonoriented_dual, node_mapper)
            # Re-validate
            t.notes = [n for n in t.notes
                       if not n.startswith("Boundaries do not isolate")
                       and "does not separate incoming/outgoing" not in n]
            is_valid, validation_notes = validate_tangle_boundaries(
                t.boundary_pairs, t.inner_nodes,
                graph, nonoriented_dual, node_mapper)
            t.boundaries_do_not_separate = not is_valid
            t.notes.extend(validation_notes)
            if not t.boundaries_do_not_separate:
                logging.info(
                    f"  Tangle {t.tangle_id}: boundaries now separate "
                    f"with passthrough pairs!")

    # Step 6: Final validation
    logging.info("Step 6: Final validation...")
    for t in tangles:
        t.notes = [n for n in t.notes
                   if not n.startswith("Boundaries do not isolate")
                   and "does not separate incoming/outgoing" not in n]
        if t.boundary_pairs and t.inner_nodes:
            is_valid, validation_notes = validate_tangle_boundaries(
                t.boundary_pairs, t.inner_nodes, graph, nonoriented_dual, node_mapper)
            t.boundaries_do_not_separate = not is_valid
            t.notes.extend(validation_notes)
        elif not t.boundary_pairs:
            t.boundaries_do_not_separate = True

        # Try to extend boundaries that fail directional separation.
        # Revert if extension doesn't fix the separation issue, so that
        # boundary changes don't mask the original no-path classification.
        if t.boundaries_do_not_separate and t.boundary_pairs:
            saved_boundary_pairs = [dict(bp) for bp in t.boundary_pairs]
            saved_inner_nodes = set(t.inner_nodes)
            saved_notes = list(t.notes)
            extend_failing_boundaries(
                t, all_scaffold_tokens, graph, nonoriented_graph, nonoriented_dual,
                cov, median_cov, node_mapper)
            if t.boundaries_do_not_separate:
                # Extension failed: revert to original boundaries.
                t.boundary_pairs = saved_boundary_pairs
                t.inner_nodes = saved_inner_nodes
                t.notes = saved_notes
                # Boundaries that don't separate mean no valid tangle can
                # be defined (per the formal definition, removing boundary
                # edges must split the component).  Classify as no_path
                # unless the tangle is already flagged multichromosomal
                # (where non-separating boundaries are expected).
                if not t.is_multichromosomal:
                    t.has_no_graph_path = True
                t.notes.append("Boundaries do not separate and extension failed")

        # Re-check graph paths (preserve no-path from extension failure)
        if not t.has_no_graph_path:
            t.notes = [n for n in t.notes if not n.startswith("No graph path")]
            for bp in t.boundary_pairs:
                path_exists, _ = find_graph_path(
                    bp['start'], bp['start_orientation'],
                    bp['end'], bp['end_orientation'],
                    graph, node_mapper)
                if not path_exists:
                    t.has_no_graph_path = True
                    t.notes.append(
                        f"No graph path from {bp['start_orientation']}{bp['start']} to "
                        f"{bp['end_orientation']}{bp['end']} "
                        f"(scaffold {bp.get('scaffold', '?')})")

    # Step 7: Merge adjacent invalid tangles that share a boundary node.
    # When two tangles from the same scaffold share a boundary and are both
    # individually invalid, merge them with outer boundaries which may form
    # a valid tangle.
    logging.info("Step 7: Merging adjacent invalid tangles with shared boundaries...")
    merged = True
    while merged:
        merged = False
        # Build index: for each scaffold, ordered list of tangle indices
        scaffold_tangle_indices = defaultdict(list)
        for i, t in enumerate(tangles):
            for g in t.gaps:
                if g.scaffold_name not in scaffold_tangle_indices or \
                        not any(idx == i for idx in [x[0] for x in scaffold_tangle_indices[g.scaffold_name]]):
                    min_gap_idx = min(gg.gap_index for gg in t.gaps if gg.scaffold_name == g.scaffold_name)
                    scaffold_tangle_indices[g.scaffold_name].append((i, min_gap_idx))
        for scaffold_name in scaffold_tangle_indices:
            scaffold_tangle_indices[scaffold_name].sort(key=lambda x: x[1])

        for scaffold_name, entries in scaffold_tangle_indices.items():
            for k in range(len(entries) - 1):
                ti_idx = entries[k][0]
                tj_idx = entries[k + 1][0]
                if ti_idx == tj_idx:
                    continue
                ti = tangles[ti_idx]
                tj = tangles[tj_idx]
                cat_i, _, is_invalid_i = classify_tangle(ti)
                cat_j, _, is_invalid_j = classify_tangle(tj)
                if not (is_invalid_i and is_invalid_j):
                    continue
                # Check if they share a boundary node
                ti_boundaries = set()
                for bp in ti.boundary_pairs:
                    ti_boundaries.add(bp['start'])
                    ti_boundaries.add(bp['end'])
                tj_boundaries = set()
                for bp in tj.boundary_pairs:
                    tj_boundaries.add(bp['start'])
                    tj_boundaries.add(bp['end'])
                shared = ti_boundaries & tj_boundaries
                if not shared:
                    continue
                # Merge: combine gaps, use outermost boundaries
                logging.info(
                    f"  Merging tangles {ti.tangle_id} and {tj.tangle_id} "
                    f"(shared boundary: {', '.join(shared)}, scaffold {scaffold_name})")
                all_cluster_gaps = ti.gaps + tj.gaps
                scaffold_gaps_map = defaultdict(list)
                for g in all_cluster_gaps:
                    scaffold_gaps_map[g.scaffold_name].append(g)
                new_boundary_pairs = []
                seen_bp = set()
                for sn, s_gaps in scaffold_gaps_map.items():
                    left_b = left_o = right_b = right_o = None
                    for g in s_gaps:
                        if g.left_boundary:
                            left_b = g.left_boundary
                            left_o = g.left_orientation
                            break
                    for g in reversed(s_gaps):
                        if g.right_boundary:
                            right_b = g.right_boundary
                            right_o = g.right_orientation
                            break
                    if left_b and right_b:
                        pkey = (left_b, right_b)
                        if pkey not in seen_bp:
                            seen_bp.add(pkey)
                            new_boundary_pairs.append({
                                'start': left_b, 'start_orientation': left_o,
                                'end': right_b, 'end_orientation': right_o,
                                'scaffold': sn,
                            })
                new_inner = ti.inner_nodes | tj.inner_nodes | shared
                if new_boundary_pairs:
                    new_inner = find_inner_nodes_from_graph(
                        new_boundary_pairs, new_inner, graph, nonoriented_dual, node_mapper)
                new_haplotypes = ti.haplotypes | tj.haplotypes
                # Check if merged is valid
                new_has_no_path = False
                new_notes = []
                for bp in new_boundary_pairs:
                    path_exists, _ = find_graph_path(
                        bp['start'], bp['start_orientation'],
                        bp['end'], bp['end_orientation'],
                        graph, node_mapper)
                    if not path_exists:
                        new_has_no_path = True
                        new_notes.append(
                            f"No graph path from {bp['start_orientation']}{bp['start']} to "
                            f"{bp['end_orientation']}{bp['end']} "
                            f"(scaffold {bp.get('scaffold', '?')})")
                bdns = False
                if new_boundary_pairs and new_inner:
                    is_valid, v_notes = validate_tangle_boundaries(
                        new_boundary_pairs, new_inner, graph, nonoriented_dual, node_mapper)
                    bdns = not is_valid
                    new_notes.extend(v_notes)
                new_tangle = DetectedTangle(
                    tangle_id=ti.tangle_id,
                    gaps=all_cluster_gaps,
                    boundary_pairs=new_boundary_pairs,
                    inner_nodes=new_inner,
                    has_no_graph_path=new_has_no_path,
                    boundaries_do_not_separate=bdns,
                    haplotypes=new_haplotypes,
                    all_scaffolds=ti.all_scaffolds | tj.all_scaffolds,
                    notes=new_notes,
                )
                new_cat, _, new_invalid = classify_tangle(new_tangle)
                if not new_invalid or new_cat != cat_i:
                    # Merge improves classification: accept
                    tangles[ti_idx] = new_tangle
                    tangles[tj_idx] = None
                    merged = True
                    logging.info(
                        f"    Accepted merge: {cat_i}+{cat_j} -> {new_cat}")
                    break
            if merged:
                break
        if merged:
            tangles = [t for t in tangles if t is not None]
            for i, t in enumerate(tangles):
                t.tangle_id = i

    logging.info(f"Final: {len(tangles)} tangles after all processing")
    return tangles


def format_tangle_report(tangles, graph, cov, node_mapper):
    """Format detected tangles as a detailed report."""
    lines = []
    lines.append(f"# Detected Tangles: {len(tangles)}")
    lines.append(f"# MIN_BOUNDARY_LENGTH={MIN_BOUNDARY_LENGTH}, "
                 f"COV_LOW_FACTOR={COV_LOW_FACTOR}, COV_HIGH_FACTOR={COV_HIGH_FACTOR}")
    lines.append("")

    for t in tangles:
        lines.append(f"=== Tangle {t.tangle_id} ===")
        lines.append(f"Haplotypes: {', '.join(sorted(t.haplotypes))}")

        if t.is_multichromosomal:
            lines.append(f"*** MULTICHROMOSOMAL ({len(t.all_scaffolds)} scaffolds) ***")
        if t.is_multihaplotype:
            lines.append(f"*** MULTIHAPLOTYPE ({len(t.haplotypes)} haplotypes) — unresolvable ***")
        if t.boundaries_not_synchronized:
            lines.append(f"*** BOUNDARIES NOT SYNCHRONIZED ***")
        if t.has_no_graph_path:
            lines.append(f"*** NO GRAPH PATH through gap ***")
        if t.has_shared_boundary:
            lines.append(f"*** SHARED BOUNDARY across haplotypes ***")
        if t.notes:
            for note in t.notes:
                lines.append(f"  NOTE: {note}")

        lines.append(f"Gaps ({len(t.gaps)}):")
        for g in t.gaps:
            lines.append(f"  scaffold: {g.scaffold_name}")
            lines.append(f"    gap_marker: {g.gap_marker}")

            if g.left_boundary:
                nid = node_mapper.get_id_for_name(g.left_boundary)
                ln = graph.nodes.get(nid, {}).get('length', '?') if nid else '?'
                cv = cov.get(nid, '?') if nid else '?'
                lines.append(f"    left_boundary:  {g.left_orientation}{g.left_boundary}  (len={ln}, cov={cv})")
            else:
                lines.append(f"    left_boundary:  NONE")

            if g.right_boundary:
                nid = node_mapper.get_id_for_name(g.right_boundary)
                ln = graph.nodes.get(nid, {}).get('length', '?') if nid else '?'
                cv = cov.get(nid, '?') if nid else '?'
                lines.append(f"    right_boundary: {g.right_orientation}{g.right_boundary}  (len={ln}, cov={cv})")
            else:
                lines.append(f"    right_boundary: NONE")

            inner_str = ', '.join(g.inner_node_names[:20])
            lines.append(f"    inner_nodes ({len(g.inner_node_names)}): {inner_str}")
            if len(g.inner_node_names) > 20:
                lines.append(f"    ... and {len(g.inner_node_names) - 20} more")

        lines.append(f"Boundary pairs ({len(t.boundary_pairs)}):")
        for bp in t.boundary_pairs:
            lines.append(f"  {bp['start_orientation']}{bp['start']}  ->  "
                         f"{bp['end_orientation']}{bp['end']}  ({bp['scaffold']})")

        lines.append(f"Total inner nodes: {len(t.inner_nodes)}")
        lines.append("")

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Detect tangles from scaffold GAF files.")
    parser.add_argument("--graph", required=True, help="Path to GFA graph file")
    parser.add_argument("--scaffolds", required=True, help="Path to scaffold GAF file")
    parser.add_argument("--coverage", required=False, help="Path to coverage CSV file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING"],
                        help="Logging level")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(levelname)s: %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(os.path.join(args.outdir, 'detect_tangles.log'), mode='w'),
        ]
    )

    node_mapper = NodeIdMapper()

    logging.info("Parsing GFA graph...")
    graph = parse_gfa(args.graph, node_mapper)
    logging.info(f"Graph: {len([n for n in graph.nodes() if n > 0])} nodes")

    if args.coverage:
        cov = read_coverage_file(args.coverage, node_mapper)
    else:
        cov = coverage_from_graph(graph)
    verify_coverage(cov, graph, node_mapper)

    median_cov = compute_median_coverage(graph, cov)
    logging.info(f"Median coverage: {median_cov:.2f}")
    logging.info(f"Boundary criteria: length >= {MIN_BOUNDARY_LENGTH}, "
                 f"coverage in [{median_cov * COV_LOW_FACTOR:.2f}, {median_cov * COV_HIGH_FACTOR:.2f}]")

    logging.info("Detecting tangles from scaffolds...")
    tangles = detect_tangles_from_scaffolds(
        args.scaffolds, graph, cov, median_cov, node_mapper)

    # Write detailed report
    report = format_tangle_report(tangles, graph, cov, node_mapper)
    report_file = os.path.join(args.outdir, 'detected_tangles.txt')
    with open(report_file, 'w') as f:
        f.write(report)
    logging.info(f"Report written to {report_file}")

    # Write JSON
    json_data = []
    for t in tangles:
        json_data.append({
            'tangle_id': t.tangle_id,
            'haplotypes': sorted(t.haplotypes),
            'is_multihaplotype': t.is_multihaplotype,
            'is_multichromosomal': t.is_multichromosomal,
            'has_no_graph_path': t.has_no_graph_path,
            'has_shared_boundary': t.has_shared_boundary,
            'boundaries_do_not_separate': t.boundaries_do_not_separate,
            'boundaries_not_synchronized': t.boundaries_not_synchronized,
            'notes': t.notes,
            'boundary_pairs': t.boundary_pairs,
            'num_gaps': len(t.gaps),
            'gaps': sorted(f"{g.scaffold_name}_gap_{g.gap_index}" for g in t.gaps),
            'all_scaffolds': sorted(t.all_scaffolds),
            'num_inner_nodes': len(t.inner_nodes),
            'inner_nodes': sorted(t.inner_nodes),
        })

    json_file = os.path.join(args.outdir, 'detected_tangles.json')
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    logging.info(f"JSON written to {json_file}")

    # Print summary to stdout
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(tangles)} tangles detected")
    print(f"{'='*60}")
    for t in tangles:
        category, flags, is_invalid = classify_tangle(t)
        pairs = "; ".join(
            f"{bp['start_orientation']}{bp['start']} -> {bp['end_orientation']}{bp['end']}"
            for bp in t.boundary_pairs)
        flag_str = f" [{', '.join(flags)}]" if flags else ""
        print(f"\nTangle {t.tangle_id}: {len(t.gaps)} gap(s), "
              f"{len(t.boundary_pairs)} boundary pair(s), "
              f"{len(t.inner_nodes)} inner nodes{flag_str}")
        gap_labels = sorted(f"{g.scaffold_name}_gap_{g.gap_index}" for g in t.gaps)
        print(f"  Gaps: {', '.join(gap_labels)}")
        if is_invalid:
            for note in t.notes:
                if note.startswith("Multichromosomal:"):
                    print(f"  >> {note}")
        else:
            print(f"  Boundaries: {pairs if pairs else 'NO VALID BOUNDARIES'}")
            if t.notes:
                for note in t.notes:
                    print(f"  >> {note}")

    # Tangle classification summary
    tangle_counts = defaultdict(int)
    for t in tangles:
        cat, _, _ = classify_tangle(t)
        tangle_counts[cat] += 1

    print(f"\n{'='*60}")
    print(f"TANGLE CLASSIFICATION SUMMARY")
    print(f"{'='*60}")
    print(f"  Valid 1-haplotype tangles:  {tangle_counts['valid_1hap']}")
    print(f"  Valid 2-haplotype tangles:  {tangle_counts['valid_2hap']}")
    print(f"  No path in graph:           {tangle_counts['no_path']}")
    print(f"  Multiscaffold (>2):         {tangle_counts['multiscaffold']}")
    print(f"  Other invalid:              {tangle_counts['other_invalid']}")
    print(f"  Total:                      {len(tangles)}")

    # Gap classification summary
    gap_counts = defaultdict(int)
    total_gaps = 0
    for t in tangles:
        cat, _, _ = classify_tangle(t)
        n_gaps = len(t.gaps)
        total_gaps += n_gaps
        gap_counts[cat] += n_gaps

    print(f"\n{'='*60}")
    print(f"GAP CLASSIFICATION SUMMARY")
    print(f"{'='*60}")
    print(f"  Valid (1-haplotype tangle): {gap_counts['valid_1hap']}")
    print(f"  Valid (2-haplotype tangle): {gap_counts['valid_2hap']}")
    print(f"  No path in graph:           {gap_counts['no_path']}")
    print(f"  Multiscaffold (>2):         {gap_counts['multiscaffold']}")
    print(f"  Other invalid:              {gap_counts['other_invalid']}")
    print(f"  Total:                      {total_gaps}")


if __name__ == '__main__':
    main()
