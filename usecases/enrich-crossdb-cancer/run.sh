#!/usr/bin/env bash
# Cross-database replication with DBRetina enrich.
#
# Seed enrich with the cancer pathways from ONE database (KEGG) and recover the
# cancer pathways that Reactome and WikiPathways curated independently, out of the
# 3090-pathway MSigDB C2:CP compendium (KEGG + Reactome + WikiPathways + BioCarta
# + PID). If cancer's shared-gene structure crosses curation boundaries, it is
# biology, not a KEGG artifact.
set -euo pipefail

DBRetina index -g data/c2cp.gmt -o c2cp --no-progress
DBRetina pairwise -i c2cp -m ochiai -c 0 --no-progress

# KEGG cancer seed + the independent (Reactome + WikiPathways) held-out truth
python build_labels.py

echo "--- pathways whose neighborhood is most enriched for the KEGG cancer seed ---"
DBRetina enrich -d c2cp_DBRetina_pairwise --module data/kegg_cancer_seed.txt \
  -m ochiai -c 20 --exclude-module -n 20 -o enrich_top.tsv
column -t -s$'\t' enrich_top.tsv | head -12

# cross-database recovery AUC + random-seed null, then figures
python validate.py
python plot.py
python extract_network.py
echo "done. see figures/"
