#!/usr/bin/env python3
"""Extract 24 cancer drugs and their target genes from DSigDB into data/drugs.gmt.

DSigDB is not shipped here. Download DSigDB_All.gmt.gz from https://dsigdb.tanlab.org/, gunzip it,
and point DSIGDB_GMT at it. data/drugs.gmt is already committed, so you can skip this and run run.sh.
For each drug we keep the DSigDB entry with the most target genes.
"""
import os

DSIGDB_GMT = os.environ.get(
    "DSIGDB_GMT",
    "/tmp/dsigdb_all.gmt",  # gunzip DSigDB_All.gmt.gz to here, or set DSIGDB_GMT
)

DRUGS = [
    "imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib",
    "gefitinib", "erlotinib", "afatinib", "lapatinib",
    "sorafenib", "sunitinib", "pazopanib", "regorafenib", "vandetanib", "cabozantinib",
    "crizotinib", "vemurafenib", "dabrafenib", "trametinib", "ruxolitinib",
    "ibrutinib", "vorinostat", "tamoxifen", "methotrexate",
]


def main():
    best = {}
    with open(DSIGDB_GMT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            name = p[0].lower()
            genes = [g.strip() for g in p[2:] if g.strip()]
            for b in DRUGS:
                if name.startswith(b) and len(genes) > len(best.get(b, ("", []))[1]):
                    best[b] = (p[0], genes)
    os.makedirs("data", exist_ok=True)
    with open("data/drugs.gmt", "w") as out:
        for b in DRUGS:
            if b in best:
                out.write(f"{b}\tna\t" + "\t".join(best[b][1]) + "\n")
    print(f"wrote data/drugs.gmt: {sum(b in best for b in DRUGS)} drugs")


if __name__ == "__main__":
    main()
