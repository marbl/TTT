# detect_tangles.py — Specification

## Purpose

Detect **tangles** in a genome assembly graph by analyzing scaffold paths. A tangle is a region of the assembly graph where scaffolding encountered unresolved complexity, marked by gap tokens (`[N...N:*]`) in scaffold GAF files. The script identifies these regions, determines their boundaries, clusters related gaps into tangles, and classifies each tangle.

## Inputs

| Argument | Required | Description |
|---|---|---|
| `--graph` | Yes | GFA assembly graph file (nodes with lengths, edges) |
| `--scaffolds` | Yes | Scaffold GAF file (tab-separated: scaffold name, path string) |
| `--coverage` | No | CSV coverage file; if omitted, coverage is derived from graph annotations |
| `--outdir` | Yes | Output directory (created if absent) |
| `--log-level` | No | `DEBUG`, `INFO` (default), or `WARNING` |

### Scaffold path format

A scaffold path is a sequence of oriented nodes (`>name`, `<name`) and gap markers (`[N<digits>N:<type>]`). Example:

```
>utig4-100<utig4-200[N50N:tangle]>utig4-300>utig4-400
```

Gap markers matching the pattern `[N\d+N:.*]` are treated as real gaps. The type label (e.g. `scaffold`, `tangle`, `rdna`) is not used for classification — all gap types are treated identically.

## Pipeline

### Step 1 — Parse scaffolds

Read the GAF file. Tokenize each scaffold path into oriented nodes and gap markers.

### Step 2 — Detect gaps and walk for boundaries

For each gap marker in each scaffold:

1. **Walk left** from the gap along the scaffold path until a node qualifies as a boundary or another gap is reached.
2. **Walk right** similarly.
3. Nodes traversed during the walk that don't qualify as boundaries are recorded as **inner nodes**.

A node qualifies as a **boundary** if:
- Its length ≥ `MIN_BOUNDARY_LENGTH` (40,000 bp)
- Its coverage is within `[median × COV_LOW_FACTOR, median × COV_HIGH_FACTOR]` (i.e., `[0.5×, 1.5×]` median)
- It exists in the graph

Coverage median is length-weighted across all graph nodes.

### Step 3 — Cluster gaps into tangles

Gaps are clustered via **union-find on shared inner nodes**: if two gaps (from any scaffolds) share an inner node, they belong to the same tangle.

**Key rule:** A gap's inner nodes are only used for clustering if a directed graph path exists between its left and right boundaries. Gaps with no graph path have walked nodes that aren't truly enclosed in a tangle region, so they must not trigger merging.

This means:
- Two gaps from the same scaffold sharing inner nodes → same tangle (e.g., adjacent gaps in utig4-908)
- Two gaps from different haplotypes sharing inner nodes → same tangle (e.g., haplotype1 and haplotype2 crossing the same graph region)
- A gap with no graph path between boundaries → stays isolated, never merges with other gaps via inner nodes

### Step 4 — Resolve shared boundaries

When a boundary node is used by gaps from different haplotypes within the same tangle, walk further outward in each scaffold to find per-haplotype distinct boundaries. The shared node becomes an inner node.

### Step 5 — Flag multichromosomal tangles

For each tangle with gaps in long scaffolds (≥ 5 Mbp):

1. Find the graph path between each boundary pair.
2. Check if any node on that path belongs to a scaffold that doesn't have a gap in this tangle.
3. If so, that scaffold **passes through** the tangle region → flag as **multichromosomal**.

### Step 6 — Final validation

For each tangle:

1. **Boundary isolation check:** Remove boundary nodes from the undirected graph. Verify that the connected component containing inner nodes is no larger than the original component minus boundaries. This confirms boundaries properly separate the tangle.
2. **Directional separation check:** For each oriented boundary node, verify it doesn't have both predecessors and successors inside the tangle component.
3. **Graph path re-check:** Confirm path existence between each boundary pair.

## Inner node computation

When boundary pairs exist, inner nodes are computed from the graph (not just the scaffold walk):

1. Remove boundary nodes from the undirected graph.
2. Find the connected component containing a seed inner node (from the walk or from BFS path).
3. All nodes in that component are inner nodes.

This captures graph nodes not mentioned in any scaffold but enclosed between boundaries.

## Tangle classification

Each tangle is classified into exactly one category:

| Category | Condition |
|---|---|
| **no_path** | No directed graph path exists between at least one boundary pair |
| **multiscaffold** | Multichromosomal (passthrough scaffolds) or multihaplotype (>2 haplotypes) |
| **other_invalid** | Boundaries don't isolate, shared boundary unresolved, or no boundaries found |
| **valid_2hap** | Valid tangle spanning exactly 2 haplotypes |
| **valid_1hap** | Valid tangle within a single haplotype |

Priority: no_path > multiscaffold > other_invalid > valid.

## Outputs

### `detected_tangles.txt`

Detailed human-readable report. For each tangle: haplotypes, flags, per-gap boundaries with lengths/coverages, inner node lists, boundary pairs.

### `detected_tangles.json`

Machine-readable JSON array. Each entry contains: tangle ID, haplotypes, all flags, boundary pairs, gap labels (`scaffold_gap_N`), scaffold list, inner node list.

### stdout summary

Per-tangle one-liner with gap count, boundary pairs, inner node count, and flags. Followed by two classification tables:

- **Tangle classification:** counts by category
- **Gap classification:** counts by category (a multi-gap tangle contributes multiple gaps to its category)

### `detect_tangles.log`

Detailed processing log written to the output directory.

## Constants

| Name | Value | Purpose |
|---|---|---|
| `MIN_BOUNDARY_LENGTH` | 40,000 | Minimum node length for boundary qualification |
| `COV_LOW_FACTOR` | 0.5 | Lower bound: `median_cov × 0.5` |
| `COV_HIGH_FACTOR` | 1.5 | Upper bound: `median_cov × 1.5` |
| `MIN_SCAFFOLD_LENGTH` | 5,000,000 | Minimum scaffold length for multichromosomal detection |

## Dependencies

- Python 3
- `networkx`
- Local modules from `src/`: `NodeIdMapper`, `parse_gfa`, `read_coverage_file`, `coverage_from_graph`, `verify_coverage`, `get_nonoriented_graph`
