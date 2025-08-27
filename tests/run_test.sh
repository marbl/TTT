#!/usr/bin/env bash
test=$1
rm tmp/*
./TTT.py --graph tests/$test/graph.noseq.gfa --coverage tests/$test/ont-coverage.csv --alignment tests/$test/alignments.gaf --outdir tmp  --boundary-nodes tests/$test/border.ids  --num-initial-paths 3
diff tmp/traversal.gaf tests/$test/expected_traversal.gaf
python3 tools/diff_compare.py tmp/traversal.gaf tests/$test/expected_traversal.gaf
