#!/usr/bin/env bash
# Type 1 vs type 2 diabetes: same name, different genetic neighbourhood.
# Run from this directory with the `dbretina` conda env active.
set -euo pipefail
PW=dm_DBRetina_pairwise

# 1. index + pairwise
DBRetina index -g data/panel.gmt -o dm --no-progress
DBRetina pairwise -i dm -m ochiai -c 0 --legacy-output --no-progress

# 2. network graph: T2D lands in the metabolic cluster, T1D in the autoimmune cluster
DBRetina graph -i dm -p $PW -m ochiai -c 0 -o dm_graph
python plot_graph.py dm_graph_edges.tsv data/labels.json figures/diabetes_graph.png 43 \
  "Type 1 and type 2 diabetes: same name, different neighbourhood (ochiai >= 43%)"

# 3. each type's nearest disorders
echo "--- type 2 diabetes' nearest ---"; DBRetina module -d $PW dis_t2d -m ochiai -c 10 --min-shared 20 -n 5
echo "--- type 1 diabetes' nearest ---"; DBRetina module -d $PW dis_t1d -m ochiai -c 10 --min-shared 20 -n 5

# 4. shared diabetes core (both types)
printf 'dis_t1d\ndis_t2d\n' > both_diabetes.txt
DBRetina geneinfo -i dm -g both_diabetes.txt -o diabetes_core

# 5. quantify the autoimmune-vs-metabolic split + figure
python analyze.py

echo "done. see figures/"
