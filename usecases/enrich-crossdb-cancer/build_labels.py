"""KEGG cancer seed + per-database / cancer labels for the C2:CP compendium.

The seed is KEGG's cancer pathways only. The held-out truth is the cancer
pathways curated independently by Reactome and WikiPathways. Both are defined
from names alone, so recovering one from the other is a fair cross-database test.
"""
import json
import re

import pyarrow.parquet as pq

names = pq.read_table("c2cp_DBRetina_pairwise/names.parquet").column("group_name").to_pylist()


def db(nm):
    m = re.match(r"(kegg|reactome|biocarta|pid|wp)_", nm)
    return m.group(1) if m else "other"


CANCER = re.compile(r"cancer|carcinoma|leukemia|melanoma|glioma")
kegg_cancer = sorted(x for x in names if db(x) == "kegg" and CANCER.search(x))
ext_cancer = sorted(x for x in names if db(x) in ("reactome", "wp") and CANCER.search(x))

with open("data/kegg_cancer_seed.txt", "w") as fh:
    fh.write("\n".join(kegg_cancer) + "\n")
with open("data/independent_cancer.txt", "w") as fh:
    fh.write("\n".join(ext_cancer) + "\n")
with open("data/labels.json", "w") as fh:
    json.dump({x: {"db": db(x), "cancer": bool(CANCER.search(x))} for x in names}, fh)

print(f"KEGG cancer seed: {len(kegg_cancer)} pathways")
print(f"independent cancer (Reactome + WikiPathways): {len(ext_cancer)} pathways")
