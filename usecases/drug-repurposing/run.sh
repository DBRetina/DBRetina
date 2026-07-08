#!/usr/bin/env bash
# Drug repurposing use case. Run from this directory with the `dbretina` conda env active.
set -euo pipefail

# 1. (optional) rebuild data/drugs.gmt from your DSigDB copy
# python build_data.py

# 2. index + pairwise
DBRetina index -g data/drugs.gmt -o drugs --no-progress
DBRetina pairwise -i drugs -m ochiai -c 0 --legacy-output --no-progress

# 3. what shares Sorafenib's targets
DBRetina neighbors -d drugs_DBRetina_pairwise "sorafenib" -m ochiai -c 0

# 4. regenerate the clustermap
mkdir -p figures && python make_figures.py

echo "done. see figures/drug_clustermap.png"
