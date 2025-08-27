#!/usr/bin/env bash
test=helo1
rm tmp/*
./TTT.py --graph tests/$test/graph.noseq.gfa --coverage tests/$test/ont-coverage.csv --boundary-nodes tests/$test/border.ids --alignment tests/$test/alignments.gaf --outdir tmp --num-initial-paths 3
diff tmp/traversal.gaf tests/$test/expected_traversal.gaf
python3 tools/diff_compare.py tmp/traversal.gaf tests/$test/expected_traversal.gaf
