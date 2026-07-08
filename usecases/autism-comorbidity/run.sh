#!/usr/bin/env bash
# Autism's coexisting conditions: comorbidity vs single-gene cause, via DBRetina.
# Run from this directory with the `dbretina` conda env active.
set -euo pipefail
PW=nz_DBRetina_pairwise

# 1. index + pairwise
DBRetina index -g data/neuro.gmt -o nz --no-progress
DBRetina pairwise -i nz -m ochiai -c 0 --legacy-output --no-progress

# 2. read the shape of the disorder network: dendrogram + force-directed graph
DBRetina export -p $PW -m ochiai --newick -o figures/nz_tree -l names
DBRetina graph -i nz -p $PW -m ochiai -c 0 -o nz_graph
python plot_graph.py nz_graph_edges.tsv data/labels.json figures/autism_graph.png 22 \
  "Autism's disease network (ochiai >= 22%)"

# 3. autism's nearest disorders (its comorbidities)
echo "--- autism's nearest disorders ---"
DBRetina module -d $PW dis_asd -m ochiai -c 10 --min-shared 20 -n 6

# 4. shared core across autism and its core comorbidities
printf 'dis_asd\ndis_epilepsy\ndis_adhd\ndis_intellectual_disability\ndis_schizophrenia\n' > core.txt
DBRetina geneinfo -i nz -g core.txt -o core_genes

# 5. rank neighbours, test polygenic-comorbidity vs control, write the figure
python analyze.py

echo "done. see figures/"
