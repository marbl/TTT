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

### Step 3b — Merge tangles with overlapping inner nodes

After initial clustering and inner node computation, tangles whose inner node sets overlap (≥80% of the larger set) are merged. The merged tangle inherits all gaps, boundary pairs, and the union of inner nodes. Inner nodes are recomputed from the graph using the dual graph after each merge. This step repeats until no more overlaps exist.

**Inflation guard:** When both tangles have >500 inner nodes (`MAX_RELIABLE_INNER`), their inner node sets may be inflated (spanning most of the graph because boundaries don't fully isolate). The merge is rejected if the combined inner node set is ≥50% of the smaller original set, indicating the overlap is an artifact of inflation rather than genuine shared interior. Tangles missed by this guard are recovered later in Step 5b via passthrough scaffold relationships.

### Step 5 — Flag multichromosomal tangles

For each tangle with gaps in long scaffolds (≥ 5 Mbp):

1. Check if any **inner node** belongs to a scaffold that doesn't have a gap in this tangle.
2. If so, that scaffold **passes through** the tangle region.
3. If the total number of scaffolds (gap + passthrough) exceeds 2, flag as **multichromosomal**.

### Step 5b — Merge related multichromosomal tangles

Multichromosomal tangles whose gap scaffolds appear as passthrough scaffolds of each other are connected by a union-find and merged into a single multichromosomal tangle. This recovers merges that Step 3b rejected due to the inflation guard — when separate tangles are actually in the same genomic region but their inflated inner node sets prevented direct merge.

### Step 5.5 — Add boundaries from passthrough scaffolds

When boundaries don't separate and a passthrough scaffold exists, the second haplotype traverses the tangle without a gap. The script discovers the passthrough scaffold's boundary pair:

1. Walk outward from the gap scaffold's start boundary until finding a node **shared** with the passthrough scaffold (appears in both scaffolds' token lists).
2. Walk outward from the end boundary similarly to find a second shared anchor.
3. Locate the shared anchors on the passthrough scaffold.
4. Between the anchors (inclusive) on the passthrough scaffold, scan for boundary-qualifying nodes.
5. Add the outermost qualifying nodes as a new boundary pair for the passthrough scaffold.
6. Recompute inner nodes and re-validate with the expanded boundary set.

The shared-anchor approach handles cases where the immediate flanking nodes of the boundary are haplotype-specific and don't exist on the passthrough scaffold. By walking further outward until a shared node is found, the algorithm reliably locates the corresponding region.

This converts 1-haplotype tangles with non-separating boundaries into valid 2-haplotype tangles when the second haplotype path exists.

### Step 6 — Final validation

For each tangle:

1. **Boundary separation check (dual graph):** Remove boundary *edges* from the nonoriented dual graph. Verify that the inner component is properly disconnected from the rest of the graph.
2. **Directional separation check:** For each oriented boundary node, verify it doesn't have both predecessors and successors inside the tangle component.
3. **Boundary extension:** If boundaries fail separation, walk further outward in the scaffold to find replacement boundaries. Revert if extension doesn't fix separation.
4. **Graph path re-check:** Confirm path existence between each boundary pair.

### Step 7 — Merge adjacent invalid tangles

When two tangles from the same scaffold share a boundary node and are both individually invalid, merge them using the outer boundaries, which may form a valid tangle. Inner nodes are recomputed via the dual graph after merging.

## Inner node computation (dual graph)

When boundary pairs exist, inner nodes are computed using the **dual graph** (not the original undirected graph):

1. Build the dual graph where vertices = junctions (canonical connections between consecutive nodes) and edges = original graph nodes.
2. Remove dual edges corresponding to boundary nodes.
3. Find the connected component containing a seed junction vertex (from a known inner node's dual edge, a BFS path node, or a boundary neighbor junction).
4. All original node IDs from edges within that component are inner nodes.

This correctly handles **node-based GFA graphs** where the tangle interior may be disconnected when boundary nodes are removed from the original graph, but remains connected through shared junctions in the dual graph. The dual graph transformation matches the paper's edge-based formalism applied to Verkko's node-based GFA format.

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

Detailed human-readable report. For each tangle: haplotypes, flags, per-gap immediate scaffold neighbors (nodes directly adjacent to the gap in the scaffold path), inner node lists, boundary pairs.

### `detected_tangles.json`

Machine-readable JSON array. Each entry contains: tangle ID, haplotypes, all flags, boundary pairs, gap labels (`scaffold_gap_N`), scaffold list, inner node list.

### Log summary

Per-tangle one-liner (via `logging.info`) with gap count, boundary pairs, inner node count, and flags. Followed by two classification tables:

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
| `MAX_RELIABLE_INNER` | 500 | Threshold above which inner node sets are considered potentially inflated |

## Dependencies

- Python 3
- `networkx`
- Local modules from `src/`: `NodeIdMapper`, `parse_gfa`, `read_coverage_file`, `coverage_from_graph`, `verify_coverage`, `get_nonoriented_graph`, `get_nonoriented_dual_graph`, `create_dual_graph`
