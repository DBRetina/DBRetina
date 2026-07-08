---
title: "How independent are the cancer Hallmark signatures?"
description: Checking whether the 50 MSigDB Hallmarks are distinct, and what the overlapping ones share.
---

The MSigDB Hallmark collection is 50 gene sets, each meant to capture one cancer-related process: apoptosis,
hypoxia, the two interferon responses, cell-cycle checkpoints, and so on. They were curated specifically to cut
the overlap and noise of older gene-set collections (Liberzon et al., 2015), and they are used across cancer
genomics as if the underlying processes (Hanahan & Weinberg, 2011) were independent. The question is how
independent they are in practice.

`index` and `pairwise` compute ochiai similarity for every pair of the 50 sets, and `export` renders the matrix as
a clustered heatmap:

```bash
DBRetina index -g h.all.gmt -o hallmark
DBRetina pairwise -i hallmark -m ochiai -c 0
DBRetina export -p hallmark_DBRetina_pairwise -m ochiai -o hallmark_heat
```

![Clustered heatmap of ochiai similarity between the 50 Hallmark gene sets](/DBRetina/img/hallmark_heatmap.png)

Most of the heatmap is dark, so the Hallmarks are mostly distinct, which is what a careful curation should
produce. A handful of pairs stand out. Ranked by ochiai:

```text
ochiai  shared  hallmark pair
 52.4     73    interferon_alpha_response  &  interferon_gamma_response
 50.5    101    estrogen_response_late     &  estrogen_response_early
 36.5     73    g2m_checkpoint             &  e2f_targets
 34.3     57    complement                 &  coagulation
 32.5     65    hypoxia                    &  glycolysis
```

## What the overlaps actually are

A high similarity is only meaningful if the shared genes are coherent. `shared-genes` lists them:

```bash
DBRetina shared-genes -d hallmark_DBRetina_pairwise -i hallmark \
    "hallmark_interferon_alpha_response" "hallmark_interferon_gamma_response"
```

The two interferon responses share 73 genes, and they are the interferon-stimulated genes: the chemokines CXCL10
and CXCL11, the antiviral effectors BST2 and EIF2AK2, the OAS family, GBP4, IFI27, and ADAR. That overlap is the
shared part of the interferon program, not a curation slip. The same check on `g2m_checkpoint` and `e2f_targets`
returns the mitotic machinery: the Aurora kinases AURKA and AURKB, CDK1 and CDK4, cyclins, the CDC25 phosphatases,
CHEK1, and survivin (BIRC5). Those two Hallmarks are two views of the cell cycle.

The Hallmarks hold up as a mostly non-overlapping set, and the few overlaps that exist are real shared biology.
`shared-genes` turns a number in the heatmap into the genes behind it.

## References

- Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database hallmark gene set collection.
  *Cell Systems* 2015;1(6):417–425. [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
- Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. *Cell* 2011;144(5):646–674.
  [doi:10.1016/j.cell.2011.02.013](https://doi.org/10.1016/j.cell.2011.02.013)
