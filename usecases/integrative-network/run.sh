#!/usr/bin/env bash
# Five-databases use case. Run from this directory with the `dbretina` conda env active.
set -euo pipefail

# 1. (optional) rebuild data/multidb.gmt from the five sources
# python build_data.py

# 2. one index over all five databases, then pairwise (containment, directional)
DBRetina index -g data/multidb.gmt -o multidb --no-progress
DBRetina pairwise -i multidb -m containment -c 0 --legacy-output --no-progress

# 3. what connects to rheumatoid arthritis: pathways (mechanism) and drugs (treatment)
echo "--- top pathways ---"
DBRetina neighbors -d multidb_DBRetina_pairwise "dis_rheumatoid_arthritis" -m containment -c 0 -n 400 \
    | grep -E "^kegg_|^reactome_|^wp_" | head -10
echo "--- top drugs ---"
DBRetina neighbors -d multidb_DBRetina_pairwise "dis_rheumatoid_arthritis" -m containment -c 0 -n 400 \
    | grep -E "^drug_" | head -10

# 4. regenerate the figure
mkdir -p figures && python make_figures.py

echo "done. see figures/multidb_ra.png"
