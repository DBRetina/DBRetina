---
title: "Alzheimer's disease and the gene-set size confound"
description: Gene-set size hides Alzheimer's true neighbours until a size-aware metric reveals the neurodegeneration signal.
---

Alzheimer's is the most-studied neurodegenerative disease, and its DisGeNET gene list is large (about 3,400 genes).
That size matters, and this case study is as much about **choosing the right metric** as about Alzheimer's itself.
The question is whether Alzheimer's shares genes with the other neurodegenerative diseases or with the vascular and
metabolic diseases so often linked to it.

The panel comprises 17 disorders drawn from DisGeNET (Piñero et al., 2017): Alzheimer's; a neurodegenerative group
(Parkinson's, ALS, Huntington's, dementia, frontotemporal dementia, Lewy body, tauopathies, prion disease); a
vascular group (ischemic stroke, cerebral infarction, vascular dementia); a metabolic group (type 2 diabetes,
obesity, metabolic syndrome); and two large controls (breast cancer, rheumatoid arthritis).

## Raw overlap tracks gene-set size

`index` and `pairwise` build the all-versus-all similarity matrix; `neighbors` ranks Alzheimer's nearest disorders
under a chosen metric.

```bash
DBRetina index -g panel.gmt -o alz
DBRetina pairwise -i alz -m ochiai -c 0
DBRetina neighbors dis_alzheimer -d alz_DBRetina_pairwise -m ochiai -c 0
```

Ranked by ochiai, the neighbour order is biologically implausible. Breast cancer, type 2 diabetes, and obesity
crowd the top, while the dementias clinically closest to Alzheimer's (tauopathies, Lewy body, frontotemporal
dementia) fall to the bottom. The ranking is tracking gene-set **size**: breast cancer's list has 6,775 genes, so
it overlaps everything, and even ochiai, which divides by the square root of the sizes, cannot fully undo the
advantage a very large set has.

![Alzheimer's neighbours ranked by ochiai versus odds ratio](/DBRetina/img/alzheimer_metric.png)

## A size-aware metric

DBRetina computes several metrics from the same comparison. The **odds ratio** measures a different quantity: how
much more overlap occurs than expected by chance, given the two sizes. Re-ranked on it (the right panel above):

```bash
DBRetina neighbors dis_alzheimer -d alz_DBRetina_pairwise -m odds_ratio -c 0
```

The order now aligns with the biology. Tauopathies lead at an odds ratio of 35, followed by vascular dementia, Lewy
body, prion disease, dementia, frontotemporal dementia, Huntington's, ALS, and Parkinson's. Breast cancer collapses
to an odds ratio of 1.0, exactly chance. Tauopathy is Alzheimer's own defining pathology, so its position at the
top is the biologically expected result under a size-aware metric.

## The odds-ratio network

The same odds-ratio scores define a network. At a cutoff of 6:

```bash
DBRetina graph -i alz -p alz_DBRetina_pairwise -m odds_ratio -c 6 -o alz_graph
```

![Force-directed network of Alzheimer's by odds ratio](/DBRetina/img/alzheimer_graph.png)

Built on the odds ratio, the network places Alzheimer's at the centre of a tight neurodegeneration and dementia
cluster, with the metabolic diseases and the two large controls at the edges.

## The shared genes

The genes Alzheimer's shares with the core dementias are the canonical ones: APP and PSEN1 (familial Alzheimer's),
MAPT (the tau gene behind tauopathies), APOE (the major risk gene), GRN (frontotemporal dementia), and SNCA (Lewy
body disease). These are the genes the size-blind ranking buried and the size-aware one restored (Bellenguez et
al., 2022).

## What this shows

The metric is not a detail. When gene sets vary widely in size, raw overlap is pulled toward the largest sets, and
a size-aware metric (odds ratio, or the p-value) is what recovers the real relationships. For Alzheimer's, that is
the difference between "looks like breast cancer" and "sits at the centre of neurodegeneration." The vascular link
is genuine as well: vascular dementia and cerebral infarction both rank high, consistent with the vascular
contribution to dementia. The broader metabolic hypothesis, however, does not rise above the neurodegeneration
signal once size is controlled.

## References

- Bellenguez C, Küçükali F, Jansen IE, et al. New insights into the genetic etiology of Alzheimer's disease and
  related dementias. *Nature Genetics* 2022;54(4):412–436.
  [doi:10.1038/s41588-022-01024-z](https://doi.org/10.1038/s41588-022-01024-z)
- Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET. *Nucleic Acids Research* 2017;45(D1):D833–D839.
  [doi:10.1093/nar/gkw943](https://doi.org/10.1093/nar/gkw943)
