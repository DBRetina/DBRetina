#!/usr/bin/env python3
"""Extract the 16-disease panel (8 autoimmune + comorbidities + controls) from DisGeNET.

DisGeNET is not shipped (see its licence). Point DISGENET_GMT at your copy and run this.
data/diseases.gmt and data/labels.json are already committed, so you can skip this and run run.sh.
"""
import json
import os

DISGENET_GMT = os.environ.get(
    "DISGENET_GMT",
    "/home/mabuelanin/dbretina-dedup/databases/DisGenet/data/disgenet.gmt",
)

# DisGeNET name -> (label used in the index, category)
DISEASES = {
    "Rheumatoid Arthritis": ("dis_ra", "autoimmune"),
    "Lupus Erythematosus, Systemic": ("dis_lupus", "autoimmune"),
    "Multiple Sclerosis": ("dis_ms", "autoimmune"),
    "Diabetes Mellitus, Insulin-Dependent": ("dis_t1d", "autoimmune"),
    "Ulcerative Colitis": ("dis_uc", "autoimmune"),
    "Crohn Disease": ("dis_crohn", "autoimmune"),
    "Psoriasis": ("dis_psoriasis", "autoimmune"),
    "Celiac Disease": ("dis_celiac", "autoimmune"),
    "Coronary Artery Disease": ("dis_coronary_artery_disease", "comorbidity"),
    "Myocardial Infarction": ("dis_myocardial_infarction", "comorbidity"),
    "Atherosclerosis": ("dis_atherosclerosis", "comorbidity"),
    "Major Depressive Disorder": ("dis_major_depression", "comorbidity"),
    "Depressive disorder": ("dis_depression", "comorbidity"),
    "Breast Carcinoma": ("dis_breast_cancer", "control"),
    "Schizophrenia": ("dis_schizophrenia", "control"),
    "Parkinson Disease": ("dis_parkinson", "control"),
}


def main():
    os.makedirs("data", exist_ok=True)
    labels = {}
    seen = set()
    with open("data/diseases.gmt", "w") as out:
        for line in open(DISGENET_GMT):
            p = line.rstrip("\n").split("\t")
            if p[0] in DISEASES and p[0] not in seen:
                label, category = DISEASES[p[0]]
                genes = [g for g in p[2:] if g.strip()]
                out.write(f"{label}\tna\t" + "\t".join(genes) + "\n")
                labels[label] = category
                seen.add(p[0])
    json.dump(labels, open("data/labels.json", "w"), indent=2)
    n_ai = sum(v == "autoimmune" for v in labels.values())
    print(f"wrote data/diseases.gmt: {len(labels)} diseases ({n_ai} autoimmune)")


if __name__ == "__main__":
    main()
