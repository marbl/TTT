#!/usr/bin/env bash
set -euo pipefail

# Run detect_tangles.py on all HPRC assemblies and summarize statistics.
#
# Usage: bash run_hprc.sh [--jobs N] [--samples SAMPLE1,SAMPLE2,...] [--dry-run]
#
# Paths:
#   Assemblies: /data/Phillippy2/projects/hprc-assemblies/assemblies-v4/<SAMPLE>/verkko-hi-c/
#   Output:     /data/antipovd2/res/TTT_paper/HPRC_stats/<SAMPLE>/

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DETECT_SCRIPT="$SCRIPT_DIR/detect_tangles.py"
ASSEMBLY_BASE="/data/Phillippy2/projects/hprc-assemblies/assemblies-v4"
OUTDIR_BASE="/data/antipovd2/res/TTT_paper/HPRC_stats"

JOBS=1
FILTER_SAMPLES=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs)    JOBS="$2"; shift 2 ;;
        --samples) FILTER_SAMPLES="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *)         echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$OUTDIR_BASE"

# Collect samples
if [[ -n "$FILTER_SAMPLES" ]]; then
    IFS=',' read -ra SAMPLES <<< "$FILTER_SAMPLES"
else
    SAMPLES=()
    for d in "$ASSEMBLY_BASE"/*/; do
        sample="$(basename "$d")"
        SAMPLES+=("$sample")
    done
fi

# Filter to those with required files
VALID_SAMPLES=()
for sample in "${SAMPLES[@]}"; do
    graph="$ASSEMBLY_BASE/$sample/verkko-hi-c/assembly.homopolymer-compressed.noseq.gfa"
    scaffolds="$ASSEMBLY_BASE/$sample/verkko-hi-c/assembly.paths.tsv"
    if [[ -f "$graph" && -f "$scaffolds" ]]; then
        VALID_SAMPLES+=("$sample")
    fi
done

echo "Found ${#VALID_SAMPLES[@]} samples with required files (out of ${#SAMPLES[@]} total)"

if $DRY_RUN; then
    echo "Dry run — would process:"
    for s in "${VALID_SAMPLES[@]}"; do echo "  $s"; done
    exit 0
fi

# Run detection
run_one() {
    local sample="$1"
    local graph="$ASSEMBLY_BASE/$sample/verkko-hi-c/assembly.homopolymer-compressed.noseq.gfa"
    local scaffolds="$ASSEMBLY_BASE/$sample/verkko-hi-c/assembly.paths.tsv"
    local outdir="$OUTDIR_BASE/$sample"

    mkdir -p "$outdir"

    if [[ -f "$outdir/detected_tangles.json" ]]; then
        echo "SKIP $sample (already done)"
        return 0
    fi

    echo "RUN  $sample ..."
    if python "$DETECT_SCRIPT" \
        --graph "$graph" \
        --scaffolds "$scaffolds" \
        --outdir "$outdir" \
        --log-level WARNING \
        > "$outdir/stdout.txt" 2> "$outdir/stderr.txt"; then
        echo "OK   $sample"
    else
        echo "FAIL $sample (exit $?)" >&2
    fi
}
export -f run_one
export ASSEMBLY_BASE OUTDIR_BASE DETECT_SCRIPT

if [[ "$JOBS" -gt 1 ]]; then
    printf '%s\n' "${VALID_SAMPLES[@]}" | \
        xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
else
    for sample in "${VALID_SAMPLES[@]}"; do
        run_one "$sample"
    done
fi

# Summarize
echo ""
echo "========================================"
echo "AGGREGATE SUMMARY"
echo "========================================"

SUMMARY_FILE="$OUTDIR_BASE/summary.tsv"
echo -e "sample\ttangles\tvalid_1hap\tvalid_2hap\tno_path\tmultiscaffold\tother_invalid\tgaps_total\tgaps_valid_1hap\tgaps_valid_2hap\tgaps_no_path\tgaps_multiscaffold\tgaps_other_invalid" > "$SUMMARY_FILE"

total_samples=0
total_tangles=0
total_valid1=0
total_valid2=0
total_nopath=0
total_multi=0
total_other=0

for sample in "${VALID_SAMPLES[@]}"; do
    stdout="$OUTDIR_BASE/$sample/stdout.txt"
    if [[ ! -f "$stdout" ]]; then
        continue
    fi

    # Parse tangle classification from stdout
    tangles=$(grep -oP 'SUMMARY: \K\d+' "$stdout" 2>/dev/null || echo 0)
    v1=$(grep 'Valid 1-haplotype tangles:' "$stdout" | head -1 | grep -oP '\d+$' || echo 0)
    v2=$(grep 'Valid 2-haplotype tangles:' "$stdout" | head -1 | grep -oP '\d+$' || echo 0)
    np=$(grep 'No path in graph:' "$stdout" | head -1 | grep -oP '\d+$' || echo 0)
    ms=$(grep 'Multiscaffold (>2):' "$stdout" | head -1 | grep -oP '\d+$' || echo 0)
    oi=$(grep 'Other invalid:' "$stdout" | head -1 | grep -oP '\d+$' || echo 0)

    # Parse gap classification
    gt=$(grep '^  Total:' "$stdout" | tail -1 | grep -oP '\d+$' || echo 0)
    gv1=$(grep 'Valid (1-haplotype tangle):' "$stdout" | grep -oP '\d+$' || echo 0)
    gv2=$(grep 'Valid (2-haplotype tangle):' "$stdout" | grep -oP '\d+$' || echo 0)
    gnp=$(grep 'No path in graph:' "$stdout" | tail -1 | grep -oP '\d+$' || echo 0)
    gms=$(grep 'Multiscaffold (>2):' "$stdout" | tail -1 | grep -oP '\d+$' || echo 0)
    goi=$(grep 'Other invalid:' "$stdout" | tail -1 | grep -oP '\d+$' || echo 0)

    echo -e "$sample\t$tangles\t$v1\t$v2\t$np\t$ms\t$oi\t$gt\t$gv1\t$gv2\t$gnp\t$gms\t$goi" >> "$SUMMARY_FILE"

    total_samples=$((total_samples + 1))
    total_tangles=$((total_tangles + tangles))
    total_valid1=$((total_valid1 + v1))
    total_valid2=$((total_valid2 + v2))
    total_nopath=$((total_nopath + np))
    total_multi=$((total_multi + ms))
    total_other=$((total_other + oi))
done

echo "Samples processed:        $total_samples"
echo "Total tangles:            $total_tangles"
echo "  Valid 1-haplotype:      $total_valid1"
echo "  Valid 2-haplotype:      $total_valid2"
echo "  No path in graph:       $total_nopath"
echo "  Multiscaffold (>2):     $total_multi"
echo "  Other invalid:          $total_other"
echo ""
echo "Per-sample summary: $SUMMARY_FILE"

# Also compute per-sample stats
echo ""
echo "Per-sample tangle count distribution:"
awk -F'\t' 'NR>1 {print $2}' "$SUMMARY_FILE" | sort -n | \
    awk '{a[NR]=$1; s+=$1} END {
        printf "  min=%d  median=%d  mean=%.1f  max=%d  (n=%d)\n",
            a[1], a[int(NR/2)+1], s/NR, a[NR], NR
    }'
