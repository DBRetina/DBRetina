#!/usr/bin/env bash
# Full DBRetina workflow: from a disease network down to the shared genes and a permutation test.
# Run from this directory with the `dbretina` conda env active.
set -euo pipefail
PW=dz_DBRetina_pairwise

# 1. index + pairwise
DBRetina index -g data/diseases.gmt -o dz --no-progress
DBRetina pairwise -i dz -m ochiai -c 0 --legacy-output --no-progress

# 2. read the shape of the network: dendrogram + force-directed graph
DBRetina export -p $PW -m ochiai --newick -o figures/dz_tree -l names
DBRetina graph -i dz -p $PW -m ochiai -c 0 -o dz_graph
python plot_graph.py dz_graph_edges.tsv data/labels.json figures/autoimmune_graph.png 40 \
  "The disease similarity network (ochiai >= 40%)"

# 3. tightest links: connected components at a high cutoff leaves only the strongest pairs
echo "--- tightest disease pairs (ochiai >= 50) ---"
DBRetina cluster -p $PW -m ochiai -c 50 -o dz50
grep -v '^#' dz50_clusters.tsv | awk -F'\t' '$2>=2'

# 4. psychiatric module via low-resolution community detection
DBRetina cluster -p $PW -m ochiai -c 30 --community --resolution 0.3 -o dzmod

# 5. rheumatoid arthritis's nearest diseases
echo "--- rheumatoid arthritis's module ---"
DBRetina module -d $PW dis_ra -m ochiai -c 20 --min-shared 30

# 6. shared core: geneinfo on the autoimmune diseases and on the controls
printf 'dis_ra\ndis_lupus\ndis_ms\ndis_t1d\ndis_uc\ndis_crohn\ndis_psoriasis\ndis_celiac\n' > autoimmune.txt
printf 'dis_breast_cancer\ndis_schizophrenia\ndis_parkinson\n' > controls.txt
DBRetina geneinfo -i dz -g autoimmune.txt -o ai_genes
DBRetina geneinfo -i dz -g controls.txt   -o ctrl_genes

# 7. prove the clustering is not a size artifact (size-matched permutation test) + figure
python prove.py

echo "done. see figures/"
