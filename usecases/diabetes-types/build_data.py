#!/usr/bin/env python3
"""Extract the diabetes panel from DisGeNET.

DisGeNET is not shipped (see its licence). Point DISGENET_GMT at your copy and run this.
data/panel.gmt and data/labels.json are already committed, so you can skip this and run run.sh.
"""
import json
import os

DISGENET_GMT = os.environ.get(
    "DISGENET_GMT",
    "/home/mabuelanin/dbretina-dedup/databases/DisGenet/data/disgenet.gmt",
)

PANEL = {
    "Diabetes Mellitus, Non-Insulin-Dependent": ("t2d", "type_2_diabetes"),
    "Diabetes Mellitus, Insulin-Dependent": ("t1d", "type_1_diabetes"),
    "Obesity": ("obesity", "metabolic"),
    "Hypertensive disease": ("hypertension", "metabolic"),
    "Metabolic Syndrome X": ("metabolic_syndrome", "metabolic"),
    "Non-alcoholic Fatty Liver Disease": ("nafld", "metabolic"),
    "Dyslipidemias": ("dyslipidemia", "metabolic"),
    "Coronary Artery Disease": ("cad", "metabolic"),
    "Atherosclerosis": ("atherosclerosis", "metabolic"),
    "Myocardial Infarction": ("mi", "metabolic"),
    "Diabetic Nephropathy": ("diabetic_nephropathy", "complication"),
    "Diabetic Retinopathy": ("diabetic_retinopathy", "complication"),
    "Rheumatoid Arthritis": ("ra", "autoimmune"),
    "Lupus Erythematosus, Systemic": ("lupus", "autoimmune"),
    "Celiac Disease": ("celiac", "autoimmune"),
    "Multiple Sclerosis": ("ms", "autoimmune"),
    "Breast Carcinoma": ("breast_cancer", "control"),
    "Schizophrenia": ("schizophrenia", "control"),
}


def main():
    os.makedirs("data", exist_ok=True)
    labels, seen = {}, set()
    with open("data/panel.gmt", "w") as out:
        for line in open(DISGENET_GMT):
            p = line.rstrip("\n").split("\t")
            if p[0] in PANEL and p[0] not in seen:
                label, category = PANEL[p[0]]
                out.write(f"dis_{label}\tna\t" + "\t".join(g for g in p[2:] if g.strip()) + "\n")
                labels[f"dis_{label}"] = category
                seen.add(p[0])
    json.dump(labels, open("data/labels.json", "w"), indent=2)
    print(f"wrote data/panel.gmt: {len(labels)} disorders")


if __name__ == "__main__":
    main()
