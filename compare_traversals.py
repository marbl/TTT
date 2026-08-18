#!/usr/bin/env python3
"""Compare TTT traversal paths with T2T reference subpaths.

For each tangle directory that has a traversal.gaf, extract the start/end
oriented nodes of each traversal, find the corresponding subpath in the T2T
reference between the same node pair, and compare.  Reverse-complement paths
are considered equal.

Usage:
    python3 compare_traversals.py \
        --ttt-runs /data/antipovd2/res/TTT_paper/giraffe/TTT_runs \
        --reference /data/antipovd2/res/TTT_paper/giraffe/T2T_utig1.paths.gaf
"""

import argparse
import os
import re
import sys
from collections import defaultdict


def parse_oriented_path(path_str):
    """Parse a path string like '>nodeA<nodeB>nodeC' into list of (orient, name)."""
    tokens = re.findall(r'([><])([^><]+)', path_str)
    return tokens  # list of (orientation, node_name)


def reverse_complement_path(path):
    """Reverse a path and flip orientations."""
    flip = {'>' : '<', '<': '>'}
    return [(flip[o], n) for o, n in reversed(path)]


def path_to_str(path):
    """Convert list of (orient, name) back to string."""
    return ''.join(f'{o}{n}' for o, n in path)


def paths_equal(p1, p2):
    """Check if two paths are equal, considering reverse complement."""
    if p1 == p2:
        return True
    if p1 == reverse_complement_path(p2):
        return True
    return False


def build_node_index(ref_paths):
    """Build index: bare_node_name -> list of (scaffold_name, position_index).

    Position index is the index in the parsed path list.
    """
    index = defaultdict(list)
    for scaffold_name, path in ref_paths.items():
        for i, (orient, name) in enumerate(path):
            index[name].append((scaffold_name, i))
    return index


def find_subpath(ref_paths, node_index, start_orient, start_name,
                 end_orient, end_name):
    """Find the subpath in reference between oriented start and end nodes.

    Returns (scaffold_name, subpath, is_reversed) or None.
    is_reversed indicates the subpath was found in reverse complement.
    """
    # Try forward: start then end in same scaffold
    candidates = []
    for scaffold_name, start_idx in node_index.get(start_name, []):
        path = ref_paths[scaffold_name]
        o_s, n_s = path[start_idx]
        if o_s != start_orient:
            continue
        # Look for end node after start
        for end_idx_offset in range(start_idx + 1, len(path)):
            o_e, n_e = path[end_idx_offset]
            if n_e == end_name and o_e == end_orient:
                subpath = path[start_idx:end_idx_offset + 1]
                candidates.append((scaffold_name, subpath, False,
                                   end_idx_offset - start_idx))
                break

    # Try reverse complement: end_rc then start_rc in same scaffold
    # If the traversal is A>B>C and the reference has <C<B<A, that's rc match
    rc_start_orient = '<' if start_orient == '>' else '>'
    rc_end_orient = '<' if end_orient == '>' else '>'
    for scaffold_name, end_idx in node_index.get(end_name, []):
        path = ref_paths[scaffold_name]
        o_e, n_e = path[end_idx]
        if o_e != rc_end_orient:
            continue
        for start_idx_offset in range(end_idx + 1, len(path)):
            o_s, n_s = path[start_idx_offset]
            if n_s == start_name and o_s == rc_start_orient:
                subpath_fwd = path[end_idx:start_idx_offset + 1]
                subpath = reverse_complement_path(subpath_fwd)
                candidates.append((scaffold_name, subpath, True,
                                   start_idx_offset - end_idx))
                break

    if not candidates:
        return None
    # Pick shortest match (closest pair)
    candidates.sort(key=lambda x: x[3])
    return candidates[0][:3]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ttt-runs',
                        default='/data/antipovd2/res/TTT_paper/giraffe/TTT_runs',
                        help='Directory containing tangle_*/manual_tangle_* subdirs')
    parser.add_argument('--subdir', default='',
                        help='Subdirectory within each tangle dir to find traversal.gaf (e.g. long_ont)')
    parser.add_argument('--reference',
                        default='/data/antipovd2/res/TTT_paper/giraffe/T2T_utig1.paths.gaf',
                        help='T2T reference paths GAF file')
    args = parser.parse_args()

    # Load reference paths
    ref_paths = {}
    with open(args.reference) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            scaffold_name = parts[0]
            path = parse_oriented_path(parts[1])
            ref_paths[scaffold_name] = path

    print(f"Loaded {len(ref_paths)} reference paths", file=sys.stderr)

    # Build node index
    node_index = build_node_index(ref_paths)

    # Process each tangle directory
    tangle_dirs = sorted([
        d for d in os.listdir(args.ttt_runs)
        if (d.startswith('tangle_') or d.startswith('manual_tangle_'))
        and os.path.isdir(os.path.join(args.ttt_runs, d))
    ])

    total = 0
    match = 0
    mismatch = 0
    not_found = 0
    t2t_gap = 0

    for tangle_dir in tangle_dirs:
        if args.subdir:
            trav_path = os.path.join(args.ttt_runs, tangle_dir, args.subdir, 'traversal.gaf')
        else:
            trav_path = os.path.join(args.ttt_runs, tangle_dir, 'traversal.gaf')
        if not os.path.exists(trav_path):
            continue

        traversals = []
        with open(trav_path) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 2:
                    continue
                trav_name = parts[0]
                trav_nodes = parse_oriented_path(parts[1])
                if len(trav_nodes) < 2:
                    continue
                traversals.append((trav_name, trav_nodes))

        for trav_name, trav_nodes in traversals:
            total += 1
            start_orient, start_name = trav_nodes[0]
            end_orient, end_name = trav_nodes[-1]

            result = find_subpath(ref_paths, node_index,
                                  start_orient, start_name,
                                  end_orient, end_name)

            if result is None:
                not_found += 1
                print(f"{tangle_dir}/{trav_name}: NO REFERENCE SUBPATH FOUND "
                      f"({start_orient}{start_name} -> {end_orient}{end_name})")
                continue

            scaffold_name, ref_subpath, is_reversed = result

            # Check if T2T subpath contains any gap
            if any('gap' in name for _, name in ref_subpath):
                t2t_gap += 1
                rc_tag = " [rc]" if is_reversed else ""
                print(f"{tangle_dir}/{trav_name}: T2T GAP "
                      f"(ref: {scaffold_name}, {len(ref_subpath)} nodes{rc_tag})")
                continue

            if paths_equal(trav_nodes, ref_subpath):
                match += 1
                rc_tag = " [rc]" if is_reversed else ""
                print(f"{tangle_dir}/{trav_name}: MATCH "
                      f"(ref: {scaffold_name}, {len(ref_subpath)} nodes{rc_tag})")
            else:
                mismatch += 1
                rc_tag = " [rc]" if is_reversed else ""
                trav_str = path_to_str(trav_nodes)
                ref_str = path_to_str(ref_subpath)
                print(f"{tangle_dir}/{trav_name}: MISMATCH "
                      f"(ref: {scaffold_name}, {len(ref_subpath)} nodes{rc_tag})")
                print(f"  TTT:  {trav_str}")
                print(f"  T2T:  {ref_str}")
                # Show node-level diff
                trav_set = set(trav_nodes)
                ref_set = set(ref_subpath)
                only_ttt = trav_set - ref_set
                only_t2t = ref_set - trav_set
                if only_ttt:
                    print(f"  Only in TTT ({len(only_ttt)}): "
                          f"{', '.join(o+n for o,n in sorted(only_ttt, key=lambda x: x[1]))}")
                if only_t2t:
                    print(f"  Only in T2T ({len(only_t2t)}): "
                          f"{', '.join(o+n for o,n in sorted(only_t2t, key=lambda x: x[1]))}")

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY: {total} traversals checked", file=sys.stderr)
    print(f"  Match:     {match}", file=sys.stderr)
    print(f"  Mismatch:  {mismatch}", file=sys.stderr)
    print(f"  T2T gap:   {t2t_gap}", file=sys.stderr)
    print(f"  Not found: {not_found}", file=sys.stderr)


if __name__ == '__main__':
    main()
