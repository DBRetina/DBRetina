---
title: "How much do independently curated pathway databases agree?"
description: KEGG, Reactome, and WikiPathways curate pathways separately; DBRetina measures how closely they agree.
---

KEGG, Reactome, and WikiPathways each build their pathways independently, from their own reading of the literature.
Whether they define a shared process with the same genes is a direct test of how consistent pathway knowledge is,
and it reduces to a comparison of gene sets.

The test set is the pathways that carry the same name across two of the databases, 87 high-confidence same-process
pairs. For each, the question is whether the matching pathway is its closest counterpart, by gene-set similarity,
in the other database.

![Cross-database pathway concordance: precision and the overlap distribution](/DBRetina/img/pathway_matching_proof.png)

The three databases agree closely. The correct same-process pathway is the closest match across databases **62%**
of the time and sits in the top five **84%** of the time. The contrast with chance is large: two pathways that name
the same process share a median **49%** of their genes, while two random pathways from different databases share a
median of **0%**. KEGG and WikiPathways agree on 98% of the cell-cycle genes and 98% of the
nucleotide-excision-repair genes, and the ErbB and chemokine signalling pathways match in the nineties.

Where the top match is not the exact counterpart, it is usually a broader or narrower version of the same process
in the other database, a genuine curation difference rather than an error.

## The similarity carries information beyond gene frequency

The similarity could in principle reflect nothing more than how common a gene is. It does not: predicting a
disease's genes from its network neighbours beats a plain gene-frequency baseline in all 16 diseases tested (mean
+0.019 AUC, Wilcoxon p < 0.001). The margin is small, but it is consistent and significant in every disease.

## What this shows

Independently curated pathway databases converge on the same genes when they describe the same biology, and
gene-set similarity captures that agreement. The same measurement locates where two databases disagree, and by how
much. The [autoimmune case study](/DBRetina/examples/autoimmune-diseases/) applies the same rigour, with a
permutation test, to a biological result.

## References

- Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 2000;28(1):27–30.
  [doi:10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway knowledgebase 2022. *Nucleic Acids Research*
  2022;50(D1):D687–D692. [doi:10.1093/nar/gkab1028](https://doi.org/10.1093/nar/gkab1028)
- Martens M, Ammar A, Riutta A, et al. WikiPathways: connecting communities. *Nucleic Acids Research*
  2021;49(D1):D613–D621. [doi:10.1093/nar/gkaa1024](https://doi.org/10.1093/nar/gkaa1024)
