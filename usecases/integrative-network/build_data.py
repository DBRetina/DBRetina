#!/usr/bin/env python3
"""Assemble data/multidb.gmt from five databases: DisGeNET + DSigDB + KEGG + Reactome + WikiPathways.

Source databases are not shipped (see their licences). Set the paths below (or the matching env vars)
and run this to rebuild. data/multidb.gmt is already committed, so you can skip this and run run.sh.
Disease names are prefixed `dis_`, drugs `drug_`; pathways keep their MSigDB source prefix (KEGG_, ...).
"""
import os

DISGENET = os.environ.get("DISGENET_GMT", "/home/mabuelanin/dbretina-dedup/databases/DisGenet/data/disgenet.gmt")
DSIGDB = os.environ.get("DSIGDB_GMT", "/tmp/dsigdb_all.gmt")
MSIGDB = os.environ.get("MSIGDB_C2", "/home/mabuelanin/dbretina-dedup/databases/MsigDB/data/c2_curated_gene_sets")

DISEASES = {
    "Rheumatoid Arthritis": "dis_rheumatoid_arthritis",
    "Lupus Erythematosus, Systemic": "dis_lupus",
    "Multiple Sclerosis": "dis_multiple_sclerosis",
    "Diabetes Mellitus, Insulin-Dependent": "dis_type1_diabetes",
    "Breast Carcinoma": "dis_breast_cancer",
}
DRUGS = ["methotrexate", "hydroxychloroquine", "sulfasalazine", "leflunomide", "tofacitinib",
         "adalimumab", "etanercept", "infliximab", "azathioprine", "prednisone", "imatinib", "sorafenib"]
PATHWAYS = ["c2.cp.kegg.v2023.1.Hs.symbols.gmt", "c2.cp.reactome.v2023.1.Hs.symbols.gmt",
            "c2.cp.wikipathways.v2023.1.Hs.symbols.gmt"]


def main():
    os.makedirs("data", exist_ok=True)
    out = open("data/multidb.gmt", "w")
    # diseases
    for line in open(DISGENET):
        p = line.rstrip("\n").split("\t")
        if p[0] in DISEASES:
            out.write(DISEASES[p[0]] + "\tna\t" + "\t".join(g for g in p[2:] if g.strip()) + "\n")
    # drugs (best DSigDB entry per drug)
    best = {}
    for line in open(DSIGDB):
        p = line.rstrip("\n").split("\t")
        nm = p[0].lower(); genes = [g.strip() for g in p[2:] if g.strip()]
        for b in DRUGS:
            if nm.startswith(b) and len(genes) > len(best.get(b, ("", []))[1]):
                best[b] = (p[0], genes)
    for b in DRUGS:
        if b in best:
            out.write("drug_" + b + "\tna\t" + "\t".join(best[b][1]) + "\n")
    # pathways (keep source prefix)
    for f in PATHWAYS:
        for line in open(os.path.join(MSIGDB, f)):
            p = line.rstrip("\n").split("\t")
            out.write(p[0].lower() + "\tna\t" + "\t".join(g for g in p[2:] if g.strip()) + "\n")
    out.close()
    print("wrote data/multidb.gmt:", sum(1 for _ in open("data/multidb.gmt")), "gene sets")


if __name__ == "__main__":
    main()
