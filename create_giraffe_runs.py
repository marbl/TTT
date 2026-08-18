#!/usr/bin/env python3
"""Generate TTT_runs subdirectories for each detected giraffe tangle."""
import os

BASE = "/data/antipovd2/res/TTT_paper/giraffe/TTT_runs"
TTT_PY = "/data/antipovd2/devel/TTT/TTT.py"
GRAPH = "../unitig-unrolled-hifi-resolved.noseq.gfa"
COVERAGE = "../unitig-unrolled-hifi-resolved.ont-coverage.csv"
ALIGNMENT = "../alns-ont.gaf"

# Tangle data from detect_tangles.py run on giraffe_detection
tangles = [
    {
        "id": 0,
        "valid": False,
        "gaps": ["haplotype1_from_utig4-1264_gap_1", "haplotype1_from_utig4-366_gap_2", "haplotype2_from_utig4-613_gap_1"],
        "reason": "MULTICHROMOSOMAL: 3 gaps from 3 different scaffolds (haplotype1_from_utig4-1264, haplotype1_from_utig4-366, haplotype2_from_utig4-613). 1073 inner nodes, 3 boundary pairs.",
    },
    {
        "id": 1,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-1264_gap_2"],
        "boundaries": [[">utig1-534", "<utig1-3407"]],
    },
    {
        "id": 2,
        "valid": False,
        "gaps": ["haplotype1_from_utig4-1367_gap_1"],
        "reason": "OTHER INVALID: 1 gap, 0 boundary pairs (no valid boundary nodes found). 5 inner nodes.",
    },
    {
        "id": 3,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-1367_gap_2"],
        "boundaries": [[">utig1-4935", "<utig1-3625"]],
    },
    {
        "id": 4,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-1367_gap_3"],
        "boundaries": [["<utig1-3625", ">utig1-3461"]],
    },
    {
        "id": 5,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-3428_gap_1", "haplotype1_from_utig4-3428_gap_2", "haplotype2_from_utig4-3427_gap_1"],
        "boundaries": [[">utig1-4018", "<utig1-1365"], ["<utig1-3237", "<utig1-3825"]],
    },
    {
        "id": 6,
        "valid": False,
        "gaps": ["haplotype1_from_utig4-4260_gap_1", "haplotype2_from_utig4-4261_gap_1"],
        "reason": "OTHER INVALID: 2 gaps from diploid scaffolds (haplotype1_from_utig4-4260, haplotype2_from_utig4-4261). 1035 inner nodes, 2 boundary pairs. Boundaries do not separate tangle from the rest of the graph.",
    },
    {
        "id": 7,
        "valid": False,
        "gaps": ["haplotype1_from_utig4-4260_gap_2"],
        "reason": "NO GRAPH PATH: 1 gap from haplotype1_from_utig4-4260. No directed graph path exists between boundary nodes. 5 inner nodes.",
    },
    {
        "id": 8,
        "valid": False,
        "gaps": [
            "haplotype1_from_utig4-759_gap_1", "haplotype1_from_utig4-759_gap_2",
            "haplotype1_from_utig4-759_gap_3", "haplotype1_from_utig4-759_gap_4",
            "haplotype1_from_utig4-759_gap_5",
            "haplotype2_from_utig4-2080_gap_1", "haplotype2_from_utig4-2080_gap_2",
            "haplotype2_from_utig4-2080_gap_3", "haplotype2_from_utig4-2080_gap_4",
            "haplotype2_from_utig4-2080_gap_5", "haplotype2_from_utig4-2080_gap_6",
            "haplotype2_from_utig4-2080_gap_7",
        ],
        "reason": "OTHER INVALID: 12 gaps from diploid scaffolds (haplotype1_from_utig4-759 x5, haplotype2_from_utig4-2080 x7). 1258 inner nodes, 2 boundary pairs. Boundaries do not separate tangle from the rest of the graph.",
    },
    {
        "id": 9,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-831_gap_2"],
        "boundaries": [[">utig1-11596", "<utig1-15334"]],
    },
    {
        "id": 10,
        "valid": False,
        "gaps": ["haplotype2_from_utig4-1815_gap_1"],
        "reason": "OTHER INVALID: 1 gap, 0 boundary pairs (no valid boundary nodes found). 7 inner nodes.",
    },
    {
        "id": 11,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-1816_gap_1", "haplotype2_from_utig4-1815_gap_2"],
        "boundaries": [["<utig1-2820", ">utig1-2971"], ["<utig1-2821", ">utig1-2974"]],
    },
    {
        "id": 12,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-759_gap_6", "haplotype2_from_utig4-2080_gap_8"],
        "boundaries": [[">utig1-12031", ">utig1-11097"], [">utig1-12030", ">utig1-11098"]],
    },
    {
        "id": 13,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-2635_gap_1", "haplotype2_from_utig4-2634_gap_1"],
        "boundaries": [[">utig1-7653", ">utig1-8446"], ["<utig1-12322", ">utig1-8447"]],
    },
    {
        "id": 14,
        "valid": True,
        "gaps": ["haplotype2_from_utig4-3427_gap_2"],
        "boundaries": [[">utig1-14647", ">utig1-2811"]],
    },
    {
        "id": 15,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-3472_gap_1", "haplotype2_from_utig4-3473_gap_1"],
        "boundaries": [["<utig1-5098", ">utig1-8565"], ["<utig1-5099", ">utig1-8566"]],
    },
    {
        "id": 16,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-366_gap_1", "haplotype2_from_utig4-367_gap_1"],
        "boundaries": [["<utig1-2854", ">utig1-2855"], ["<utig1-2853", ">utig1-2856"]],
    },
    {
        "id": 17,
        "valid": False,
        "gaps": ["haplotype2_from_utig4-3844_gap_1", "haplotype2_from_utig4-3844_gap_3"],
        "reason": "MULTICHROMOSOMAL: 2 gaps from haplotype2_from_utig4-3844 (non-adjacent, gap_1 and gap_3). 121 inner nodes, 1 boundary pair. Multiple scaffolds pass through tangle boundaries.",
    },
    {
        "id": 18,
        "valid": True,
        "gaps": ["haplotype2_from_utig4-3844_gap_2"],
        "boundaries": [["<utig1-1676", ">utig1-1574"]],
    },
    {
        "id": 19,
        "valid": True,
        "gaps": ["haplotype2_from_utig4-3844_gap_4"],
        "boundaries": [[">utig1-3390", ">utig1-4056"]],
    },
    {
        "id": 20,
        "valid": True,
        "gaps": [
            "haplotype1_from_utig4-4100_gap_1", "haplotype1_from_utig4-4100_gap_2",
            "haplotype2_from_utig4-4099_gap_1", "haplotype2_from_utig4-4099_gap_2",
        ],
        "boundaries": [["<utig1-6745", "<utig1-15026"], ["<utig1-4599", "<utig1-15027"]],
    },
    {
        "id": 21,
        "valid": True,
        "gaps": ["haplotype1_from_utig4-831_gap_1", "haplotype2_from_utig4-830_gap_1"],
        "boundaries": [["<utig1-2434", ">utig1-7447"], ["<utig1-2435", ">utig1-7446"]],
    },
    {
        "id": 22,
        "valid": True,
        "gaps": ["na_unused_utig4-1514_gap_1"],
        "boundaries": [[">utig1-2344", "<utig1-1965"]],
    },
    {
        "id": 23,
        "valid": True,
        "gaps": ["na_unused_utig4-1514_gap_2", "na_unused_utig4-3760_gap_1"],
        "boundaries": [["<utig1-9636", "<utig1-8236"], ["<utig1-9637", "<utig1-8235"]],
    },
]

for t in tangles:
    dirname = f"tangle_{t['id']:02d}"
    dirpath = os.path.join(BASE, dirname)
    os.makedirs(dirpath, exist_ok=True)

    if t["valid"]:
        # Write boundary_nodes.txt
        bn_path = os.path.join(dirpath, "boundary_nodes.txt")
        with open(bn_path, "w") as f:
            for pair in t["boundaries"]:
                f.write(f"{pair[0]} {pair[1]}\n")

        # Write run.sh
        run_path = os.path.join(dirpath, "run.sh")
        with open(run_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Tangle {t['id']}: {', '.join(t['gaps'])}\n")
            f.write(f"# Boundaries: {'; '.join(p[0] + ' -> ' + p[1] for p in t['boundaries'])}\n\n")
            f.write(f"python3 {TTT_PY} \\\n")
            f.write(f"  --graph {GRAPH} \\\n")
            f.write(f"  --coverage {COVERAGE} \\\n")
            f.write(f"  --alignment {ALIGNMENT} \\\n")
            f.write(f"  --boundary-nodes boundary_nodes.txt \\\n")
            f.write(f"  --outdir .\n")
        os.chmod(run_path, 0o755)
    else:
        # Write info.txt with explanation
        info_path = os.path.join(dirpath, "info.txt")
        with open(info_path, "w") as f:
            f.write(f"Tangle {t['id']} - INVALID\n")
            f.write(f"Gaps: {', '.join(t['gaps'])}\n")
            f.write(f"Reason: {t['reason']}\n")

    print(f"Created {dirpath} ({'valid' if t['valid'] else 'INVALID'})")

print(f"\nDone: {sum(1 for t in tangles if t['valid'])} valid, {sum(1 for t in tangles if not t['valid'])} invalid")
