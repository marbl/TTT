#!/usr/bin/env python3
"""
Standalone tangle detection from scaffold GAF files.
Parses scaffolds for [N...N:gap] markers, finds boundary nodes by walking
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
from src.input_parsing import parse_gfa, read_coverage_file, coverage_from_graph, verify_coverage

# Constants
MIN_BOUNDARY_LENGTH = 50_000       # Minimum node length to qualify as boundary
COV_LOW_FACTOR      = 0.5          # Boundary cov >= median * COV_LOW_FACTOR
COV_HIGH_FACTOR     = 1.5          # Boundary cov <= median * COV_HIGH_FACTOR
GAP_PATTERN         = re.compile(r'^\[N\d+N:.*\]$')


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
    is_multihaplotype: bool = False   # >2 distinct haplotypes
    is_multichromosomal: bool = False  # >2 scaffolds traverse this region
    has_no_graph_path: bool = False    # no path through gap in graph
    has_shared_boundary: bool = False  # same boundary in multiple haplotypes
    boundaries_do_not_separate: bool = False  # boundaries don't isolate tangle from graph
    boundaries_not_synchronized: bool = False  # left/right boundaries not from same graph fork
    haplotypes: set = field(default_factory=set)
    all_scaffolds: set = field(default_factory=set)  # all scaffolds touching this region
    notes: list = field(default_factory=list)


def tokenize_scaffold_path(path_str):
    """Parse a scaffold path string into tokens handling <> nodes and [...] gaps."""
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


def is_real_gap(token):
    """Check if token is a [N...N:gap] marker (NOT gapont)."""
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
    if length < MIN_BOUNDARY_LENGTH:
        return False
    if coverage < median_cov * COV_LOW_FACTOR:
        return False
    if coverage > median_cov * COV_HIGH_FACTOR:
        return False
    return True


def is_in_graph(node_name, node_mapper, graph):
    """Check if node_name exists in graph."""
    if not node_mapper.has_name(node_name):
        return False
    node_id = node_mapper.get_id_for_name(node_name)
    return node_id in graph.nodes or (-node_id) in graph.nodes


def get_haplotype(scaffold_name):
    """Extract haplotype identifier from scaffold name (e.g. 'haplotype1' from 'haplotype1_from_utig4-961')."""
    # Common patterns: haplotype1_from_..., haplotype2_from_..., na_from_...
    parts = scaffold_name.split('_from_')
    return parts[0] if parts else scaffold_name


def oriented_node_id(node_name, orientation, node_mapper):
    """Get the signed node id from name + orientation (<  means negative, > means positive)."""
    if not node_mapper.has_name(node_name):
        return None
    nid = node_mapper.get_id_for_name(node_name)
    if orientation == '<':
        return -nid
    return nid


def check_graph_path_exists(left_name, left_orient, right_name, right_orient,
                            graph, node_mapper, max_depth=200):
    """
    Check if there is a directed path in the graph from the left boundary node
    to the right boundary node. Uses BFS with limited depth.
    """
    start_id = oriented_node_id(left_name, left_orient, node_mapper)
    end_id = oriented_node_id(right_name, right_orient, node_mapper)
    if start_id is None or end_id is None:
        return False
    if start_id not in graph.nodes or end_id not in graph.nodes:
        return False

    visited = set()
    queue = deque()
    for succ in graph.successors(start_id):
        queue.append((succ, 1))
        visited.add(succ)

    while queue:
        node, depth = queue.popleft()
        if node == end_id:
            return True
        if depth >= max_depth:
            continue
        for succ in graph.successors(node):
            if succ not in visited:
                visited.add(succ)
                queue.append((succ, depth + 1))

    return False


def validate_tangle_boundaries(boundary_pairs, inner_node_names, graph, node_mapper):
    """
    Validate that boundary nodes properly separate the tangle from the rest of the graph.
    Mirrors the logic in new_identify_tangle_nodes from input_parsing.py:
    1. Build undirected graph, remove ALL boundary nodes together
    2. Check inner nodes form an isolated connected component
    3. Check no boundary node has both predecessors AND successors inside the
       isolated tangle component in its scaffold orientation

    Returns (is_valid, notes) where:
      - is_valid: True if boundaries properly separate
      - notes: list of issue descriptions
    """
    notes = []

    # Collect all boundary node IDs (positive)
    boundary_node_ids = set()
    boundary_orientations = {}  # node_id -> set of oriented_ids from scaffold
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

    # Collect all inner node IDs (positive)
    inner_node_ids = set()
    for name in inner_node_names:
        if node_mapper.has_name(name):
            inner_node_ids.add(node_mapper.get_id_for_name(name))

    if not boundary_node_ids or not inner_node_ids:
        return True, notes  # nothing to validate

    # Build undirected graph (positive IDs only, like get_nonoriented_graph)
    indirect_graph = nx.MultiGraph()
    for node in graph.nodes():
        if node < 0:
            continue
        indirect_graph.add_node(node, **graph.nodes[node])
    for u, v in graph.edges():
        indirect_graph.add_edge(abs(u), abs(v))

    # Find an inner node that's in the graph
    inside_node = None
    for nid in inner_node_ids:
        if nid in indirect_graph:
            inside_node = nid
            break
    if inside_node is None:
        notes.append("No inner nodes found in graph")
        return False, notes

    original_component = nx.node_connected_component(indirect_graph, inside_node)

    # Remove ALL boundary nodes at once
    for bnid in boundary_node_ids:
        if bnid in indirect_graph:
            indirect_graph.remove_node(bnid)

    if inside_node not in indirect_graph:
        notes.append("Inner node was removed (was it also a boundary?)")
        return False, notes

    tangle_component = nx.node_connected_component(indirect_graph, inside_node)

    # Check: the tangle component must be smaller than original
    # (all boundaries together must actually isolate the tangle)
    boundaries_isolate = (len(tangle_component) + len(boundary_node_ids) < len(original_component))
    is_valid = True
    if not boundaries_isolate:
        notes.append("Boundaries do not isolate tangle from the rest of the connected component")
        is_valid = False

    # Check no boundary node has both predecessors AND successors inside
    # the tangle component in its oriented direction from the scaffold.
    # If oriented_id has successors in tangle, it's incoming.
    # If oriented_id has predecessors in tangle, it's outgoing.
    # Having BOTH means it doesn't cleanly separate.
    for bnid, oriented_ids in boundary_orientations.items():
        for oriented_id in oriented_ids:
            if oriented_id not in graph.nodes:
                continue
            has_pred_in = any(abs(p) in tangle_component for p in graph.predecessors(oriented_id))
            has_succ_in = any(abs(s) in tangle_component for s in graph.successors(oriented_id))
            if has_pred_in and has_succ_in:
                bname = node_mapper.node_id_to_name_safe(bnid)
                orient_char = '>' if oriented_id > 0 else '<'
                notes.append(
                    f"Boundary {orient_char}{bname} does not separate incoming/outgoing")
                is_valid = False

    return is_valid, notes


def walk_for_boundary(tokens, gap_idx, direction, graph, cov, median_cov, node_mapper):
    """
    Walk from gap_idx in given direction (+1=right, -1=left).
    Returns (boundary_node_name, orientation, inner_nodes) or (None, '', inner_nodes).
    Stops if another gap marker is hit.
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
    for gi in gap_indices:
        gap = GapInfo(scaffold_name=scaffold_name, gap_marker=tokens[gi].raw)

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


def cluster_gaps_into_tangles(all_gaps, graph, node_mapper):
    """Cluster gaps into tangles using union-find on shared nodes.
    Also detects: multihaplotype tangles, shared boundaries, no-path gaps."""
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
        all_names = set(gap.inner_node_names)
        if gap.left_boundary:
            all_names.add(gap.left_boundary)
        if gap.right_boundary:
            all_names.add(gap.right_boundary)
        for name in all_names:
            node_to_gaps[name].append(i)

    for name, gap_indices in node_to_gaps.items():
        for j in range(1, len(gap_indices)):
            union(gap_indices[0], gap_indices[j])

    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)

    tangles = []
    for tid, (root, indices) in enumerate(sorted(clusters.items())):
        cluster_gaps = [all_gaps[i] for i in indices]

        inner = set()
        for g in cluster_gaps:
            inner.update(g.inner_node_names)

        # Collect haplotypes
        haplotypes = set()
        for g in cluster_gaps:
            haplotypes.add(get_haplotype(g.scaffold_name))

        notes = []

        # Check for >2 haplotypes
        is_multihaplotype = len(haplotypes) > 2
        if is_multihaplotype:
            notes.append(f"Multihaplotype ({len(haplotypes)} haplotypes: {', '.join(sorted(haplotypes))}): unresolvable")

        # Check for shared boundaries: same node used as boundary in different haplotypes
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

        # Build boundary pairs (deduplicated)
        seen_pairs = set()
        boundary_pairs = []
        for g in cluster_gaps:
            if g.left_boundary and g.right_boundary:
                pair_key = (g.left_boundary, g.right_boundary)
                if pair_key not in seen_pairs:
                    seen_pairs.add(pair_key)
                    boundary_pairs.append({
                        'start': g.left_boundary,
                        'start_orientation': g.left_orientation,
                        'end': g.right_boundary,
                        'end_orientation': g.right_orientation,
                        'scaffold': g.scaffold_name,
                    })

        # Validate boundaries separate the tangle from graph
        # (mirrors new_identify_tangle_nodes logic)
        boundaries_do_not_separate = False
        if boundary_pairs:
            is_valid, validation_notes = validate_tangle_boundaries(
                boundary_pairs, inner, graph, node_mapper)
            if not is_valid:
                boundaries_do_not_separate = True
            notes.extend(validation_notes)

        # Check graph path existence for each boundary pair
        has_no_graph_path = False
        for bp in boundary_pairs:
            path_exists = check_graph_path_exists(
                bp['start'], bp['start_orientation'],
                bp['end'], bp['end_orientation'],
                graph, node_mapper)
            if not path_exists:
                has_no_graph_path = True
                notes.append(
                    f"No graph path from {bp['start_orientation']}{bp['start']} to "
                    f"{bp['end_orientation']}{bp['end']} "
                    f"(scaffold {bp.get('scaffold', '?')})")
        # Also flag gaps with no boundaries at all
        for g in cluster_gaps:
            if g.left_boundary is None and g.right_boundary is None:
                has_no_graph_path = True
                notes.append(f"Both boundaries missing (scaffold {g.scaffold_name})")

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


def find_passthrough_boundaries(tangle_inner_nodes, tangle_boundary_names,
                                 tangle_gap_scaffolds,
                                 all_scaffold_tokens, graph, cov, median_cov, node_mapper):
    """
    For a tangle detected from gaps, find other scaffolds that pass through
    the same inner nodes WITHOUT a gap. These represent the other haplotype(s)
    going straight through the tangle region.

    Returns list of dicts: {'start', 'start_orientation', 'end', 'end_orientation', 'scaffold'}
    """
    inner_set = set(tangle_inner_nodes)
    if not inner_set:
        return []

    passthrough_pairs = []

    for scaffold_name, tokens in all_scaffold_tokens.items():
        # Skip scaffolds that already contributed gaps to this tangle
        if scaffold_name in tangle_gap_scaffolds:
            continue

        # Find token indices that overlap with inner nodes
        overlap_indices = []
        for i, t in enumerate(tokens):
            if not t.is_gap and t.node_name in inner_set:
                overlap_indices.append(i)

        if not overlap_indices:
            continue

        # Check there's no gap within this overlap region
        min_idx = min(overlap_indices)
        max_idx = max(overlap_indices)
        has_gap_in_region = any(
            is_real_gap(tokens[i]) for i in range(min_idx, max_idx + 1))
        if has_gap_in_region:
            continue

        # Walk left from min_idx to find left boundary
        left_boundary = None
        left_orient = ''
        idx = min_idx - 1
        while idx >= 0:
            t = tokens[idx]
            if is_real_gap(t):
                break
            if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                if node_qualifies_as_boundary(t.node_name, graph, cov, median_cov, node_mapper):
                    left_boundary = t.node_name
                    left_orient = t.orientation
                    break
            idx -= 1

        # Walk right from max_idx to find right boundary
        right_boundary = None
        right_orient = ''
        idx = max_idx + 1
        while idx < len(tokens):
            t = tokens[idx]
            if is_real_gap(t):
                break
            if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                if node_qualifies_as_boundary(t.node_name, graph, cov, median_cov, node_mapper):
                    right_boundary = t.node_name
                    right_orient = t.orientation
                    break
            idx += 1

        if left_boundary and right_boundary:
            passthrough_pairs.append({
                'start': left_boundary,
                'start_orientation': left_orient,
                'end': right_boundary,
                'end_orientation': right_orient,
                'scaffold': scaffold_name,
            })
            logging.debug(
                f"  Passthrough in {scaffold_name}: "
                f"{left_orient}{left_boundary} -> {right_orient}{right_boundary}")

    return passthrough_pairs


def resolve_shared_boundaries(tangle, all_scaffold_tokens, graph, cov, median_cov, node_mapper):
    """
    When a boundary node is shared across haplotypes (e.g. both haplotypes use
    <utig1-25358 as left boundary), walk one step further away from the tangle
    in each scaffold to find per-haplotype distinct boundaries.

    The shared node becomes an inner node, and the new (per-scaffold) nodes
    before/after it become the boundaries.
    """
    if not tangle.has_shared_boundary:
        return

    # Find which boundary nodes are shared
    boundary_to_haplotypes = defaultdict(set)
    boundary_to_pairs = defaultdict(list)  # node_name -> list of (pair_idx, side)
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
        # The shared node will become inner
        tangle.inner_nodes.add(shared_node)

        # For each boundary pair that uses this shared node, look in its scaffold
        # to find the next qualifying boundary node further away from tangle
        for pair_idx, side in boundary_to_pairs[shared_node]:
            bp = tangle.boundary_pairs[pair_idx]
            scaffold_name = bp['scaffold']
            if scaffold_name not in all_scaffold_tokens:
                continue

            tokens = all_scaffold_tokens[scaffold_name]

            # Find where the shared node appears in this scaffold
            shared_token_indices = [
                i for i, t in enumerate(tokens)
                if not t.is_gap and t.node_name == shared_node
            ]
            if not shared_token_indices:
                continue

            # Determine walk direction: for 'start' (left boundary), walk left (-1)
            # For 'end' (right boundary), walk right (+1)
            if side == 'start':
                # Walk further left (away from tangle)
                direction = -1
                shared_idx = shared_token_indices[0]
            else:
                # Walk further right (away from tangle)
                direction = +1
                shared_idx = shared_token_indices[-1]

            # Walk from shared_idx in the direction away from tangle
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


def _get_undirected_neighbors(node_name, graph, node_mapper):
    """Get all undirected neighbors of a node (predecessors + successors of both orientations)."""
    if not node_mapper.has_name(node_name):
        return set()
    nid = node_mapper.get_id_for_name(node_name)
    neighbors = set()
    for s in graph.successors(nid):
        neighbors.add(abs(s))
    for p in graph.predecessors(nid):
        neighbors.add(abs(p))
    for s in graph.successors(-nid):
        neighbors.add(abs(s))
    for p in graph.predecessors(-nid):
        neighbors.add(abs(p))
    neighbors.discard(nid)  # don't include self
    return neighbors


def _are_synchronized(name1, name2, graph, node_mapper):
    """Check if two boundary nodes are synchronized (share a common graph neighbor or are neighbors)."""
    if name1 == name2:
        return True
    n1 = _get_undirected_neighbors(name1, graph, node_mapper)
    n2 = _get_undirected_neighbors(name2, graph, node_mapper)
    nid1 = node_mapper.get_id_for_name(name1) if node_mapper.has_name(name1) else None
    nid2 = node_mapper.get_id_for_name(name2) if node_mapper.has_name(name2) else None
    return bool(n1 & n2) or ((nid1 is not None and nid2 is not None) and (nid2 in n1))


def fix_flipped_boundary_pairs(tangle, graph, node_mapper):
    """
    Detect and fix boundary pairs where two haplotypes traverse the tangle in
    opposite directions (one pair is the reverse of the other).

    For example:
      BP1: >A -> >B (hap2, left-to-right)
      BP2: <B'-> <A' (hap1, right-to-left)
    If start1 syncs with end2 AND end1 syncs with start2, BP2 is flipped.
    Fix by swapping BP2's start/end so both pairs go in the same direction.
    """
    if len(tangle.boundary_pairs) != 2:
        return False
    if len(tangle.haplotypes) != 2:
        return False

    bp1, bp2 = tangle.boundary_pairs[0], tangle.boundary_pairs[1]
    if get_haplotype(bp1['scaffold']) == get_haplotype(bp2['scaffold']):
        return False

    # Check if already synchronized normally
    start_synced = _are_synchronized(bp1['start'], bp2['start'], graph, node_mapper)
    end_synced = _are_synchronized(bp1['end'], bp2['end'], graph, node_mapper)
    if start_synced and end_synced:
        return False  # already fine

    # Check cross-sync: start1 with end2, end1 with start2
    cross_start = _are_synchronized(bp1['start'], bp2['end'], graph, node_mapper)
    cross_end = _are_synchronized(bp1['end'], bp2['start'], graph, node_mapper)

    if cross_start and cross_end:
        logging.info(
            f"  Tangle {tangle.tangle_id}: fixing flipped BP2 "
            f"({bp2['scaffold']}): swapping start/end")
        # Swap BP2's start and end
        old_start = bp2['start']
        old_start_orient = bp2['start_orientation']
        bp2['start'] = bp2['end']
        bp2['start_orientation'] = bp2['end_orientation']
        bp2['end'] = old_start
        bp2['end_orientation'] = old_start_orient
        return True

    return False


def fix_unsynchronized_boundaries(tangle, all_scaffold_tokens, graph, cov,
                                   median_cov, node_mapper, max_walk_steps=5):
    """
    For 2-haplotype tangles where boundaries are not synchronized, try to fix
    by walking further outward in each scaffold to find boundaries that ARE
    synchronized with the other haplotype's boundary.

    For each unsynchronized side (start or end):
    - Take the pair that is closer to the tangle (shorter boundary node) and
      walk further outward in its scaffold
    - For each candidate boundary node encountered, check if it synchronizes
      with the other haplotype's boundary
    - If found, update the boundary pair and move old boundary to inner nodes
    """
    if len(tangle.boundary_pairs) != 2:
        return False
    if len(tangle.haplotypes) != 2:
        return False

    bp1, bp2 = tangle.boundary_pairs[0], tangle.boundary_pairs[1]
    if get_haplotype(bp1['scaffold']) == get_haplotype(bp2['scaffold']):
        return False

    changed = False

    for side, orient_key, direction in [
        ('start', 'start_orientation', -1),  # walk left for start boundary
        ('end', 'end_orientation', +1),       # walk right for end boundary
    ]:
        name1, name2 = bp1[side], bp2[side]
        if _are_synchronized(name1, name2, graph, node_mapper):
            continue  # already synced

        # Try walking further in EACH scaffold to find a synced boundary
        for bp_to_fix, bp_other in [(bp1, bp2), (bp2, bp1)]:
            other_name = bp_other[side]
            scaffold_name = bp_to_fix['scaffold']
            if scaffold_name not in all_scaffold_tokens:
                continue

            tokens = all_scaffold_tokens[scaffold_name]
            current_name = bp_to_fix[side]

            # Find current boundary position in scaffold
            boundary_indices = [
                i for i, t in enumerate(tokens)
                if not t.is_gap and t.node_name == current_name
            ]
            if not boundary_indices:
                continue

            start_idx = boundary_indices[0] if direction == -1 else boundary_indices[-1]

            # Walk further outward
            idx = start_idx + direction
            steps = 0
            while 0 <= idx < len(tokens) and steps < max_walk_steps:
                t = tokens[idx]
                if is_real_gap(t):
                    break
                if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                    if node_qualifies_as_boundary(t.node_name, graph, cov,
                                                   median_cov, node_mapper):
                        # Check if this candidate syncs with the other haplotype
                        if _are_synchronized(t.node_name, other_name, graph, node_mapper):
                            logging.info(
                                f"  Tangle {tangle.tangle_id}: walked {side} "
                                f"boundary of {scaffold_name} from "
                                f"{bp_to_fix[orient_key]}{current_name} to "
                                f"{t.orientation}{t.node_name} (syncs with "
                                f"{other_name})")
                            # Old boundary becomes inner
                            tangle.inner_nodes.add(current_name)
                            bp_to_fix[side] = t.node_name
                            bp_to_fix[orient_key] = t.orientation
                            changed = True
                            break
                        else:
                            steps += 1
                    else:
                        # Non-qualifying node — add to inner and continue
                        tangle.inner_nodes.add(t.node_name)
                idx += direction

            if changed:
                break  # Fixed this side, move on

    return changed


def check_boundaries_synchronized(tangle, graph, node_mapper):
    """
    For 2-haplotype tangles with exactly 2 boundary pairs, check that the
    left (start) boundaries share a common graph neighbor, and separately
    that the right (end) boundaries share a common graph neighbor.

    Two boundary nodes are "synchronized" if they are siblings — they share
    a predecessor or successor in the undirected graph (i.e. they fork from
    or merge into a common node).

    For example: left boundaries utig1-47927, utig1-47928 are synchronized
    if both have a common predecessor like utig1-47926. But utig1-47927 and
    utig1-52436 are NOT synchronized (they come from different parts of graph).
    """
    if len(tangle.boundary_pairs) != 2:
        return
    if len(tangle.haplotypes) != 2:
        return

    bp1, bp2 = tangle.boundary_pairs[0], tangle.boundary_pairs[1]
    # Only check if from different haplotypes
    if get_haplotype(bp1['scaffold']) == get_haplotype(bp2['scaffold']):
        return

    # Check start (left) boundaries
    start1, start2 = bp1['start'], bp2['start']
    if not _are_synchronized(start1, start2, graph, node_mapper):
        tangle.boundaries_not_synchronized = True
        tangle.notes.append(
            f"Left boundaries not synchronized: {start1} and {start2} "
            f"share no common graph neighbor")

    # Check end (right) boundaries
    end1, end2 = bp1['end'], bp2['end']
    if not _are_synchronized(end1, end2, graph, node_mapper):
        tangle.boundaries_not_synchronized = True
        tangle.notes.append(
            f"Right boundaries not synchronized: {end1} and {end2} "
            f"share no common graph neighbor")


def _test_isolation(boundary_names, inner_node_ids, graph, node_mapper):
    """Quick test: do given boundary names isolate the inner nodes?
    Returns (isolates, component_size)."""
    boundary_ids = set()
    for name in boundary_names:
        if node_mapper.has_name(name):
            boundary_ids.add(node_mapper.get_id_for_name(name))

    if not boundary_ids or not inner_node_ids:
        return False, 0

    # Build undirected graph
    ug = nx.Graph()
    for node in graph.nodes():
        if node > 0:
            ug.add_node(node)
    for u, v in graph.edges():
        ug.add_edge(abs(u), abs(v))

    # Find an inner node in the graph
    inside_node = None
    for nid in inner_node_ids:
        if nid in ug:
            inside_node = nid
            break
    if inside_node is None:
        return False, 0

    orig_size = len(nx.node_connected_component(ug, inside_node))

    # Remove boundaries
    for bid in boundary_ids:
        if bid in ug:
            ug.remove_node(bid)

    if inside_node not in ug:
        return False, 0

    comp = nx.node_connected_component(ug, inside_node)
    isolates = (len(comp) + len(boundary_ids) < orig_size)
    return isolates, len(comp)


def fix_non_isolating_boundaries(tangle, all_scaffold_tokens, graph, cov,
                                  median_cov, node_mapper, max_walk_steps=10):
    """
    For tangles whose boundaries don't isolate the inner nodes from the rest
    of the graph, try walking further outward in each scaffold to find
    boundaries that DO isolate.

    For each boundary pair, for each side (end first, then start):
    - Test if current boundaries isolate
    - If not, walk outward in the scaffold past the current boundary
    - For each qualifying candidate, tentatively replace the boundary and
      test isolation
    - If isolation achieved, keep the new boundary
    - If walk exhausted without finding isolation, REVERT to original boundary
    """
    if not tangle.boundaries_do_not_separate:
        return False
    if not tangle.boundary_pairs:
        return False

    inner_ids = set()
    for name in tangle.inner_nodes:
        if node_mapper.has_name(name):
            inner_ids.add(node_mapper.get_id_for_name(name))

    def get_all_boundary_names():
        names = set()
        for bp in tangle.boundary_pairs:
            names.add(bp['start'])
            names.add(bp['end'])
        return names

    changed = False

    for bp in tangle.boundary_pairs:
        scaffold_name = bp['scaffold']
        if scaffold_name not in all_scaffold_tokens:
            continue

        tokens = all_scaffold_tokens[scaffold_name]

        # Try end side first (direction=+1), then start side (direction=-1)
        for side, orient_key, direction in [
            ('end', 'end_orientation', +1),
            ('start', 'start_orientation', -1),
        ]:
            # Test if current boundaries already isolate
            all_bnames = get_all_boundary_names()
            current_inner_ids = set(inner_ids)
            for name in tangle.inner_nodes:
                if node_mapper.has_name(name):
                    current_inner_ids.add(node_mapper.get_id_for_name(name))

            isolates, _ = _test_isolation(all_bnames, current_inner_ids, graph, node_mapper)
            if isolates:
                continue

            # Save original state for revert
            orig_name = bp[side]
            orig_orient = bp[orient_key]
            added_inner = []

            current_name = bp[side]

            boundary_indices = [
                i for i, t in enumerate(tokens)
                if not t.is_gap and t.node_name == current_name
            ]
            if not boundary_indices:
                continue

            start_idx = boundary_indices[0] if direction == -1 else boundary_indices[-1]

            idx = start_idx + direction
            steps = 0
            found = False
            while 0 <= idx < len(tokens) and steps < max_walk_steps:
                t = tokens[idx]
                if is_real_gap(t):
                    break
                if not t.is_gap and is_in_graph(t.node_name, node_mapper, graph):
                    if node_qualifies_as_boundary(t.node_name, graph, cov,
                                                   median_cov, node_mapper):
                        # Tentatively set new boundary
                        old_name = bp[side]
                        bp[side] = t.node_name
                        bp[orient_key] = t.orientation
                        added_inner.append(old_name)
                        tangle.inner_nodes.add(old_name)
                        if node_mapper.has_name(old_name):
                            current_inner_ids.add(node_mapper.get_id_for_name(old_name))

                        new_bnames = get_all_boundary_names()
                        new_isolates, new_comp = _test_isolation(
                            new_bnames, current_inner_ids, graph, node_mapper)

                        if new_isolates:
                            logging.info(
                                f"  Tangle {tangle.tangle_id}: walked {side} "
                                f"boundary of {scaffold_name} from "
                                f"{orig_name} to {t.node_name} "
                                f"(now isolates, comp={new_comp})")
                            changed = True
                            found = True
                            break
                        else:
                            steps += 1
                    else:
                        added_inner.append(t.node_name)
                        tangle.inner_nodes.add(t.node_name)
                        if node_mapper.has_name(t.node_name):
                            current_inner_ids.add(node_mapper.get_id_for_name(t.node_name))
                idx += direction

            if not found:
                # Revert boundary to original
                bp[side] = orig_name
                bp[orient_key] = orig_orient
                # Remove added inner nodes (they were speculative)
                for n in added_inner:
                    tangle.inner_nodes.discard(n)

    return changed


def detect_tangles_from_scaffolds(scaffold_file, graph, cov, median_cov, node_mapper):
    """Full pipeline: parse scaffolds, find gaps, walk boundaries,
    find passthrough haplotypes, cluster."""

    # Step 1: Parse and tokenize ALL scaffolds
    all_scaffold_tokens = {}  # scaffold_name -> tokens
    all_scaffold_paths = {}   # scaffold_name -> path_str
    with open(scaffold_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            scaffold_name = parts[0]
            path_str = parts[1]
            all_scaffold_tokens[scaffold_name] = tokenize_scaffold_path(path_str)
            all_scaffold_paths[scaffold_name] = path_str

    # Step 2: Detect gaps in each scaffold
    all_gaps = []
    for scaffold_name, path_str in all_scaffold_paths.items():
        gaps = detect_gaps_in_scaffold(
            scaffold_name, path_str, graph, cov, median_cov, node_mapper)
        all_gaps.extend(gaps)

    logging.info(f"Found {len(all_gaps)} gaps across all scaffolds")
    logging.info (f"Total scaffolds processed: {len(all_scaffold_paths)}")
    # Step 3: Initial clustering
    tangles = cluster_gaps_into_tangles(all_gaps, graph, node_mapper)
    logging.info(f"Clustered into {len(tangles)} tangles (before passthrough search)")

    # Step 4: For each tangle, find other scaffolds passing through without gap
    for t in tangles:
        tangle_gap_scaffolds = set(g.scaffold_name for g in t.gaps)
        tangle_boundary_names = set()
        for bp in t.boundary_pairs:
            tangle_boundary_names.add(bp['start'])
            tangle_boundary_names.add(bp['end'])

        passthrough_pairs = find_passthrough_boundaries(
            t.inner_nodes, tangle_boundary_names, tangle_gap_scaffolds,
            all_scaffold_tokens, graph, cov, median_cov, node_mapper)

        if passthrough_pairs:
            for pp in passthrough_pairs:
                hap = get_haplotype(pp['scaffold'])
                t.haplotypes.add(hap)
                logging.info(
                    f"  Tangle {t.tangle_id}: passthrough from {hap} "
                    f"({pp['scaffold']}): "
                    f"{pp['start_orientation']}{pp['start']} -> "
                    f"{pp['end_orientation']}{pp['end']}")

            # Merge passthrough boundary pairs (deduplicated)
            existing_keys = set(
                (bp['start'], bp['end']) for bp in t.boundary_pairs)
            for pp in passthrough_pairs:
                key = (pp['start'], pp['end'])
                if key not in existing_keys:
                    existing_keys.add(key)
                    t.boundary_pairs.append(pp)

            # Re-validate boundaries with the expanded set
            t.notes = [n for n in t.notes
                       if not n.startswith("Boundaries do not isolate")
                       and "does not separate incoming/outgoing" not in n]
            if t.boundary_pairs:
                is_valid, validation_notes = validate_tangle_boundaries(
                    t.boundary_pairs, t.inner_nodes, graph, node_mapper)
                t.boundaries_do_not_separate = not is_valid
                t.notes.extend(validation_notes)

            # Re-check shared boundaries
            boundary_to_haplotypes = defaultdict(set)
            for bp in t.boundary_pairs:
                hap = get_haplotype(bp['scaffold'])
                boundary_to_haplotypes[bp['start']].add(hap)
                boundary_to_haplotypes[bp['end']].add(hap)
            shared = [nd for nd, haps in boundary_to_haplotypes.items() if len(haps) > 1]
            t.has_shared_boundary = bool(shared)
            # Remove old shared-boundary notes and add fresh ones
            t.notes = [n for n in t.notes if not n.startswith("Shared boundary")]
            for sn in shared:
                haps = sorted(boundary_to_haplotypes[sn])
                t.notes.append(f"Shared boundary node {sn} used by haplotypes: {', '.join(haps)}")

            # Re-check multihaplotype
            t.is_multihaplotype = len(t.haplotypes) > 2
            if t.is_multihaplotype:
                t.notes = [n for n in t.notes if not n.startswith("Multihaplotype")]
                t.notes.append(
                    f"Multihaplotype ({len(t.haplotypes)} haplotypes: "
                    f"{', '.join(sorted(t.haplotypes))}): unresolvable")

            # Re-check graph paths for new pairs
            t.has_no_graph_path = False
            t.notes = [n for n in t.notes if not n.startswith("No graph path")]
            for bp in t.boundary_pairs:
                path_exists = check_graph_path_exists(
                    bp['start'], bp['start_orientation'],
                    bp['end'], bp['end_orientation'],
                    graph, node_mapper)
                if not path_exists:
                    t.has_no_graph_path = True
                    t.notes.append(
                        f"No graph path from {bp['start_orientation']}{bp['start']} to "
                        f"{bp['end_orientation']}{bp['end']} "
                        f"(scaffold {bp.get('scaffold', '?')})")
            for g in t.gaps:
                if g.left_boundary is None and g.right_boundary is None:
                    t.has_no_graph_path = True

    # Step 5: Resolve shared boundaries by walking further in scaffolds
    logging.info("Step 5: Resolving shared boundaries...")
    for t in tangles:
        if t.has_shared_boundary:
            resolve_shared_boundaries(t, all_scaffold_tokens, graph, cov, median_cov, node_mapper)

    # Step 6: Collect all scaffolds touching each tangle (gaps + passthroughs)
    # and flag multichromosomal tangles (>2 scaffolds)
    logging.info("Step 6: Checking for multichromosomal tangles...")
    for t in tangles:
        t.all_scaffolds = set(g.scaffold_name for g in t.gaps)
        for bp in t.boundary_pairs:
            t.all_scaffolds.add(bp['scaffold'])
        if len(t.all_scaffolds) > 2:
            t.is_multichromosomal = True
            t.notes.append(
                f"Multichromosomal: {len(t.all_scaffolds)} scaffolds traverse this region: "
                f"{', '.join(sorted(t.all_scaffolds))}")
            logging.info(
                f"  Tangle {t.tangle_id}: multichromosomal "
                f"({len(t.all_scaffolds)} scaffolds)")

    # Step 7: Fix boundary synchronization for 2-haplotype tangles
    # 7a: Detect and fix flipped boundary pairs
    logging.info("Step 7a: Detecting flipped boundary pairs...")
    for t in tangles:
        if not t.is_multichromosomal:
            fix_flipped_boundary_pairs(t, graph, node_mapper)

    # 7b: Walk further in scaffolds to find synchronized boundaries
    logging.info("Step 7b: Walking further to synchronize boundaries...")
    for t in tangles:
        if not t.is_multichromosomal:
            fix_unsynchronized_boundaries(
                t, all_scaffold_tokens, graph, cov, median_cov, node_mapper)

    # 7c: Final boundary sync check (flag remaining unsynchronized)
    logging.info("Step 7c: Checking boundary synchronization...")
    for t in tangles:
        check_boundaries_synchronized(t, graph, node_mapper)

    # Step 7d: Fix non-isolating boundaries by walking further
    logging.info("Step 7d: Fixing non-isolating boundaries...")
    for t in tangles:
        if t.boundaries_do_not_separate:
            fix_non_isolating_boundaries(
                t, all_scaffold_tokens, graph, cov, median_cov, node_mapper)

    # Step 8: Final re-validation of all boundaries and paths
    logging.info("Step 8: Final validation...")
    for t in tangles:
        # Re-validate boundary separation
        t.notes = [n for n in t.notes
                   if not n.startswith("Boundaries do not isolate")
                   and "does not separate incoming/outgoing" not in n]
        if t.boundary_pairs and t.inner_nodes:
            is_valid, validation_notes = validate_tangle_boundaries(
                t.boundary_pairs, t.inner_nodes, graph, node_mapper)
            t.boundaries_do_not_separate = not is_valid
            t.notes.extend(validation_notes)
        elif not t.boundary_pairs:
            t.boundaries_do_not_separate = True

        # Re-check graph paths
        t.has_no_graph_path = False
        t.notes = [n for n in t.notes if not n.startswith("No graph path")]
        for bp in t.boundary_pairs:
            path_exists = check_graph_path_exists(
                bp['start'], bp['start_orientation'],
                bp['end'], bp['end_orientation'],
                graph, node_mapper)
            if not path_exists:
                t.has_no_graph_path = True
                t.notes.append(
                    f"No graph path from {bp['start_orientation']}{bp['start']} to "
                    f"{bp['end_orientation']}{bp['end']} "
                    f"(scaffold {bp.get('scaffold', '?')})")
        for g in t.gaps:
            if g.left_boundary is None and g.right_boundary is None:
                t.has_no_graph_path = True

    logging.info(f"Final: {len(tangles)} tangles after all processing")
    return tangles


def format_tangle_report(tangles, graph, cov, node_mapper):
    """Format detected tangles as a readable report."""
    lines = []
    lines.append(f"# Detected Tangles: {len(tangles)}")
    lines.append(f"# MIN_BOUNDARY_LENGTH={MIN_BOUNDARY_LENGTH}, "
                 f"COV_LOW_FACTOR={COV_LOW_FACTOR}, COV_HIGH_FACTOR={COV_HIGH_FACTOR}")
    lines.append("")

    for t in tangles:
        lines.append(f"=== Tangle {t.tangle_id} ===")
        lines.append(f"Haplotypes: {', '.join(sorted(t.haplotypes))}")

        # Show flags
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

    # Write human-readable report
    report = format_tangle_report(tangles, graph, cov, node_mapper)
    report_file = os.path.join(args.outdir, 'detected_tangles.txt')
    with open(report_file, 'w') as f:
        f.write(report)
    logging.info(f"Report written to {report_file}")

    # Write JSON
    json_data = []
    for t in tangles:
        tdata = {
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
            'scaffolds': sorted(set(g.scaffold_name for g in t.gaps)),
            'all_scaffolds': sorted(t.all_scaffolds),
            'num_inner_nodes': len(t.inner_nodes),
            'inner_nodes': sorted(t.inner_nodes),
        }
        json_data.append(tdata)

    json_file = os.path.join(args.outdir, 'detected_tangles.json')
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    logging.info(f"JSON written to {json_file}")

    # Print summary to stdout
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(tangles)} tangles detected")
    print(f"{'='*60}")
    for t in tangles:
        scaff = sorted(t.all_scaffolds)
        pairs = "; ".join(
            f"{bp['start_orientation']}{bp['start']} -> {bp['end_orientation']}{bp['end']}"
            for bp in t.boundary_pairs)
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
        flag_str = f" [{', '.join(flags)}]" if flags else ""
        print(f"\nTangle {t.tangle_id}: {len(t.gaps)} gap(s), "
              f"{len(t.boundary_pairs)} boundary pair(s), "
              f"{len(t.inner_nodes)} inner nodes{flag_str}")
        print(f"  Haplotypes: {', '.join(sorted(t.haplotypes))}")
        print(f"  Scaffolds: {', '.join(scaff)}")
        print(f"  Boundaries: {pairs if pairs else 'NO VALID BOUNDARIES'}")
        if t.notes:
            for note in t.notes:
                print(f"  >> {note}")


if __name__ == '__main__':
    #while debugging, expected run  python detect_tangles.py --graph HG002_detection/assembly.homopolymer-compressed.noseq.gfa --scaffolds HG002_detection/scaff_rukki.paths.gaf --outdir rerun_tst
    main()
