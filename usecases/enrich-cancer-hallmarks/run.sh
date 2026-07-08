#!/usr/bin/env bash
# Recovering the Hallmarks of Cancer from KEGG by guilt-by-association.
#
# Told only the names of specific cancers (breast, glioma, melanoma, the
# leukaemias, ...), DBRetina enrich ranks every other pathway by how enriched
# its gene-sharing neighborhood is for that module. The pathways it surfaces are
# the oncogenic signaling network -- the Hallmarks of Cancer.
set -euo pipefail

DBRetina index -g data/kegg.gmt -o kegg --no-progress
DBRetina pairwise -i kegg -m ochiai -c 0 --no-progress

# objective module: the KEGG cancer-type pathways, defined by NAME only
python build_module.py

echo "--- pathways whose neighborhood is most enriched for the cancer module ---"
DBRetina enrich -d kegg_DBRetina_pairwise --module data/cancer_module.txt \
  -m ochiai -c 20 --exclude-module -n 20 -o enrich_top.tsv
column -t -s$'\t' enrich_top.tsv | head -16

# held-out recovery (leave-one-out) + a size-matched random-module null, then figures
python validate.py
python plot.py
# the interactive tutorial network (nodes/edges/holdout ranks) -> network.json
python extract_network.py
echo "done. see figures/ and network.json"
