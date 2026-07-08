---
title: "One name, two mechanisms: type 1 and type 2 diabetes"
description: Type 1 diabetes groups with autoimmune disease and type 2 with metabolic disease, from gene lists alone.
---

Type 1 and type 2 diabetes share a name and a symptom, high blood glucose, but not a mechanism. Type 1 is an
autoimmune disease, in which the immune system destroys the insulin-producing beta cells. Type 2 is metabolic, in
which tissues stop responding to insulin, usually alongside obesity. If gene-set overlap tracks biology rather than
shared nomenclature, the two types should fall into different company despite the common name.

The panel comprises 18 disorders drawn from DisGeNET (Piñero et al., 2017): the two diabetes types; a metabolic
group (obesity, hypertension, metabolic syndrome, fatty liver, dyslipidaemia, coronary artery disease,
atherosclerosis, myocardial infarction); the two diabetic complications (nephropathy, retinopathy); an autoimmune
group (rheumatoid arthritis, lupus, celiac, multiple sclerosis); and two controls.

## Building the network

`index` and `pairwise` build the similarity matrix, and `graph` renders it as a node-and-edge network. At an ochiai
cutoff of 43:

```bash
DBRetina index -g panel.gmt -o dm
DBRetina pairwise -i dm -m ochiai -c 0
DBRetina graph -i dm -p dm_DBRetina_pairwise -m ochiai -c 43 -o dm_graph
```

![Force-directed network of the diabetes panel](/DBRetina/img/diabetes_graph.png)

The separation is clear. Type 2 diabetes sits inside the metabolic cluster, connected to obesity, hypertension,
metabolic syndrome, and the cardiovascular diseases. Type 1 diabetes sits across the graph inside the autoimmune
cluster, next to multiple sclerosis, lupus, celiac, and rheumatoid arthritis. The two share a name but fall into
opposite neighbourhoods.

## Each type's nearest disorders

`module` lists each type's characteristic neighbours as a set:

```bash
DBRetina module -d dm_DBRetina_pairwise dis_t2d -m ochiai -c 10 --min-shared 20 -n 5
DBRetina module -d dm_DBRetina_pairwise dis_t1d -m ochiai -c 10 --min-shared 20 -n 5
```

Type 2 diabetes returns obesity, hypertension, atherosclerosis, metabolic syndrome, and myocardial infarction, all
metabolic. Type 1 diabetes returns rheumatoid arthritis, multiple sclerosis, and lupus, all autoimmune.

## Quantifying the split

Similarity to each group puts a number on the split.

![Type 1 leans autoimmune, type 2 leans metabolic](/DBRetina/img/diabetes_split.png)

Type 1 diabetes is more similar to the autoimmune diseases than to the metabolic ones (42% versus 35%), and type 2
is the mirror image (43% metabolic versus 33% autoimmune). The word "diabetes" points in two different genetic
directions.

## What they share, and what they don't

`geneinfo` on both diabetes types together returns the genes they hold in common:

```bash
DBRetina geneinfo -i dm -g both_diabetes.txt -o diabetes_core
```

The two types share 969 genes, the diabetes core of glucose and insulin handling: INS, TCF7L2, GCK, KCNJ11, HNF1A,
WFS1, SLC30A8. That shared core is why they carry the same name. What separates them is the company each keeps.
Type 1 adds the autoimmune risk loci it shares with the other autoimmune diseases (HLA-DRB1, HLA-DQB1, CTLA4,
IL2RA, PTPN22), and type 2 adds the metabolic and obesity genes it shares with the metabolic diseases (PPARG,
ADIPOQ, LEP, LEPR, FTO). These are the two gene families that distinguish the diseases in the clinic (Barrett et
al., 2009; Mahajan et al., 2018).

## What this shows

From gene lists alone, two diseases that share a name separate into their mechanistic families: autoimmune for type
1, metabolic for type 2. The shared glucose core and the divergent risk genes both fall out of the same pairwise
comparison, with no clinical annotation supplied.

## References

- Barrett JC, Clayton DG, Concannon P, et al. Genome-wide association study and meta-analysis find that over 40
  loci affect risk of type 1 diabetes. *Nature Genetics* 2009;41(6):703–707.
  [doi:10.1038/ng.381](https://doi.org/10.1038/ng.381)
- Mahajan A, Taliun D, Thurner M, et al. Fine-mapping type 2 diabetes loci to single-variant resolution using
  high-density imputation and islet-specific epigenome maps. *Nature Genetics* 2018;50(11):1505–1513.
  [doi:10.1038/s41588-018-0241-6](https://doi.org/10.1038/s41588-018-0241-6)
- Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET. *Nucleic Acids Research* 2017;45(D1):D833–D839.
  [doi:10.1093/nar/gkw943](https://doi.org/10.1093/nar/gkw943)
