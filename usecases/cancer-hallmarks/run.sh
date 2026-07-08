#!/usr/bin/env bash
# Cancer Hallmarks use case. Run from this directory with the `dbretina` conda env active.
# Needs the MSigDB Hallmark GMT (download from https://www.gsea-msigdb.org/gsea/msigdb/).
set -euo pipefail

HALLMARK_GMT="${HALLMARK_GMT:-/home/mabuelanin/dbretina-dedup/databases/MsigDB/data/h_hallmark_gene_sets/h.all.v2023.1.Hs.symbols.gmt}"

# 1. index + pairwise
DBRetina index -g "$HALLMARK_GMT" -o hallmark --no-progress
DBRetina pairwise -i hallmark -m ochiai -c 0 --no-progress

# 2. the heatmap (DBRetina's export draws it)
DBRetina export -p hallmark_DBRetina_pairwise -m ochiai -o hallmark_heat
mkdir -p figures && cp hallmark_heat_heatmap.png figures/hallmark_heatmap.png

# 3. what the overlapping Hallmark pairs share
echo "--- interferon alpha vs gamma ---"
DBRetina shared-genes -d hallmark_DBRetina_pairwise -i hallmark \
    "hallmark_interferon_alpha_response" "hallmark_interferon_gamma_response"
echo "--- g2m checkpoint vs e2f targets ---"
DBRetina shared-genes -d hallmark_DBRetina_pairwise -i hallmark \
    "hallmark_g2m_checkpoint" "hallmark_e2f_targets"

echo "done. see figures/hallmark_heatmap.png"
