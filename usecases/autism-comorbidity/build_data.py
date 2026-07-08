#!/usr/bin/env python3
"""Extract the 18-disorder autism panel from DisGeNET.

DisGeNET is not shipped (see its licence). Point DISGENET_GMT at your copy and run this.
data/neuro.gmt and data/labels.json are already committed, so you can skip this and run run.sh.
"""
import json
import os

DISGENET_GMT = os.environ.get(
    "DISGENET_GMT",
    "/home/mabuelanin/dbretina-dedup/databases/DisGenet/data/disgenet.gmt",
)

# DisGeNET name -> (label used in the index, category)
PANEL = {
    "Autism Spectrum Disorders": ("asd", "autism"),
    "Intellectual Disability": ("intellectual_disability", "neurodev"),
    "Epilepsy": ("epilepsy", "neurodev"),
    "Attention deficit hyperactivity disorder": ("adhd", "neurodev"),
    "Global developmental delay": ("dev_delay", "neurodev"),
    "Schizophrenia": ("schizophrenia", "psychiatric"),
    "Bipolar Disorder": ("bipolar", "psychiatric"),
    "Major Depressive Disorder": ("depression", "psychiatric"),
    "Anxiety Disorders": ("anxiety", "psychiatric"),
    "Obsessive-Compulsive Disorder": ("ocd", "psychiatric"),
    "Gilles de la Tourette syndrome": ("tourette", "psychiatric"),
    "Rett Syndrome": ("rett", "syndromic"),
    "Fragile X Syndrome": ("fragile_x", "syndromic"),
    "Tuberous Sclerosis": ("tuberous_sclerosis", "syndromic"),
    "Down Syndrome": ("down_syndrome", "syndromic"),
    "Rheumatoid Arthritis": ("ra", "control"),
    "Breast Carcinoma": ("breast_cancer", "control"),
    "Coronary Artery Disease": ("cad", "control"),
}


def main():
    os.makedirs("data", exist_ok=True)
    labels = {}
    seen = set()
    with open("data/neuro.gmt", "w") as out:
        for line in open(DISGENET_GMT):
            p = line.rstrip("\n").split("\t")
            if p[0] in PANEL and p[0] not in seen:
                label, category = PANEL[p[0]]
                genes = [g for g in p[2:] if g.strip()]
                out.write(f"dis_{label}\tna\t" + "\t".join(genes) + "\n")
                labels[f"dis_{label}"] = category
                seen.add(p[0])
    json.dump(labels, open("data/labels.json", "w"), indent=2)
    print(f"wrote data/neuro.gmt: {len(labels)} disorders")


if __name__ == "__main__":
    main()
