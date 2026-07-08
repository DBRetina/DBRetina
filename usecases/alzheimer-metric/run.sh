#!/usr/bin/env bash
# Alzheimer's: why the metric matters. A size-aware metric recovers the real neurodegeneration signal.
# Run from this directory with the `dbretina` conda env active.
set -euo pipefail
PW=alz_DBRetina_pairwise

# 1. index + pairwise (all metrics are computed together)
DBRetina index -g data/panel.gmt -o alz --no-progress
DBRetina pairwise -i alz -m ochiai -c 0 --legacy-output --no-progress

# 2. the trap: raw overlap ranks the biggest sets first
echo "--- Alzheimer's neighbours by ochiai (size-biased) ---"
DBRetina neighbors dis_alzheimer -d $PW -m ochiai -c 0 -n 6

# 3. the fix: a size-aware metric recovers the neurodegeneration signal
echo "--- Alzheimer's neighbours by odds ratio (size-aware) ---"
DBRetina neighbors dis_alzheimer -d $PW -m odds_ratio -c 0 -n 6

# 4. network on the odds ratio: Alzheimer's inside the neurodegeneration cluster
DBRetina graph -i alz -p $PW -m odds_ratio -c 0 -o alz_graph
python plot_graph.py alz_graph_edges.tsv data/labels.json figures/alzheimer_graph.png 6 \
  "Alzheimer's network by odds ratio (OR >= 6)" inv

# 5. before/after figure
python analyze.py

echo "done. see figures/"
