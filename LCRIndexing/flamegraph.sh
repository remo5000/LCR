#!/bin/sh
# - Assumes that perf.data exists
# - Outputs a svg in this folder

FLAME_GRAPH_DIR=~/FlameGraph
#PERF_SCRIPT_OUT=$(mktemp)

echo "perf script";
perf script > out.perf;
echo "flamegraph out";
cat out.perf \
	| $FLAME_GRAPH_DIR/stackcollapse-perf.pl \
	| $FLAME_GRAPH_DIR/flamegraph.pl \
	> out.svg;
echo "flamegraph done";

