"""Define the objective cancer module + class labels from KEGG pathway names.

The module is defined from names ALONE (no gene information), so recovering it
from gene overlap is a fair test, not a circular one.
"""
import json
import re

import pyarrow.parquet as pq

names = pq.read_table("kegg_DBRetina_pairwise/names.parquet").column("group_name").to_pylist()
cancer = sorted(x for x in names if re.search(r"cancer|carcinoma|leukemia|melanoma|glioma", x))

with open("data/cancer_module.txt", "w") as fh:
    fh.write("\n".join(cancer) + "\n")

labels = {x: ("cancer" if x in set(cancer) else "other") for x in names}
with open("data/labels.json", "w") as fh:
    json.dump(labels, fh, indent=0)

print(f"cancer module: {len(cancer)} pathways -> data/cancer_module.txt")
