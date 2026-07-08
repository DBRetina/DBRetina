---
title: "The shared genetic core of autoimmune disease"
description: A full DBRetina workflow, from a disease network down to the shared genes and a permutation test.
---

Autoimmune diseases target different tissues and present as unrelated conditions in the clinic. Rheumatoid
arthritis attacks the joints, type 1 diabetes the pancreas, multiple sclerosis the nerves, and Crohn's disease the
gut. They nonetheless share genetic risk, much of it in immune and HLA genes, and that sharing has been mapped
across the autoimmune diseases by large genetic studies (Cotsapas et al., 2011). The analysis that follows tests
whether the same structure is recoverable from disease-gene lists alone, and runs end to end, from the disease
network to the autoimmune-specific genes and a permutation test against a size-matched null.

The panel comprises 16 diseases drawn from DisGeNET (Piñero et al., 2017): eight autoimmune, and eight
non-autoimmune diseases (cardiovascular, psychiatric, cancer, Parkinson's) as controls. A genuine signal should be
specific to the autoimmune group rather than a side effect of large gene lists, which is what the controls are
there to test.

## Building the network

`index` builds the hashed representation of the gene sets, and `pairwise` computes all-versus-all ochiai similarity
with the cutoff at zero so that no edges are dropped.

```bash
DBRetina index -g diseases.gmt -o dz
DBRetina pairwise -i dz -m ochiai -c 0
```

## Reading the structure

`export` renders the pairwise matrix as a dendrogram, which shows the structure before any cutoff is imposed.

```bash
DBRetina export -p dz_DBRetina_pairwise -m ochiai --newick -o dz_tree -l names
```

![Circular dendrogram of the 16 diseases](/DBRetina/img/disease_dendrogram.png)

The dendrogram is computed from gene-set overlap alone, with no disease labels supplied, and it recovers the
expected neighbourhoods. The eight autoimmune diseases form one arc, with inflammatory bowel disease (Crohn's and
ulcerative colitis) and the rheumatoid arthritis and lupus pair sitting closest. The cardiovascular diseases
(myocardial infarction, atherosclerosis, coronary artery disease) form another, the psychiatric diseases
(depression and schizophrenia) a third, and breast cancer and Parkinson's sit apart from the rest.

The `graph` command represents the same relationships as a node-and-edge network. Retaining edges above 40%
similarity exposes the two dense neighbourhoods and the bridge between them:

```bash
DBRetina graph -i dz -p dz_DBRetina_pairwise -m ochiai -c 40 -o dz_graph
```

![Force-directed network of the 16 diseases](/DBRetina/img/autoimmune_graph.png)

The autoimmune diseases form one connected cluster, with the thick Crohn's and ulcerative colitis link standing
out; the cardiovascular diseases form a tight triangle of their own (atherosclerosis, myocardial infarction, and
coronary artery disease); and rheumatoid arthritis sits between the two as the busiest hub.

## The tightest links

The dendrogram shows structure; `cluster` quantifies it. In connected-components mode, raising the cutoff removes
all but the strongest edges, so the tightest disease pairs are the ones that survive last.

```bash
DBRetina cluster -p dz_DBRetina_pairwise -m ochiai -c 50 -o dz50
```

At ochiai 50, the two surviving clusters are the two textbook pairs: `{crohn, uc}`, the inflammatory bowel
diseases, and `{ra, lupus}`, the systemic autoimmune diseases. At a cutoff of 55 only the inflammatory bowel pair
remains, the single strongest relationship in the network. Community detection at low resolution instead separates
the psychiatric module, `{depression, schizophrenia}`, as its own community, apart from the somatic diseases.

## Rheumatoid arthritis's neighbours

`module` lists a single disease's characteristic neighbours as a set. For rheumatoid arthritis:

```bash
DBRetina module -d dz_DBRetina_pairwise dis_ra -m ochiai -c 20 --min-shared 30
```

The nearest diseases are led by lupus and Crohn's disease, both autoimmune. The output is a set of disease names,
which feeds directly into the next step.

## The shared gene core

`geneinfo` reports, for every gene, how many groups in a list carry it. Run once on the autoimmune diseases and
once on the controls, the difference defines the autoimmune-specific core: genes shared widely across autoimmunity
and absent from the control diseases.

```bash
DBRetina geneinfo -i dz -g autoimmune.txt -o ai_genes
DBRetina geneinfo -i dz -g controls.txt   -o ctrl_genes
```

The genes shared by seven or eight of the eight autoimmune diseases and carried by none of the controls are
UBASH3A, IFIH1, CLEC16A, CCR1, BACH2, and ANKRD55, followed by TAGAP, PRKCQ, and MASP2. Each is a known
autoimmune-risk gene, and their absence from the controls is what makes the core specific rather than a set of
genes common to disease lists in general.

## Ruling out a size artifact

Large gene lists overlap to some degree simply because they are large, so the clustering requires a null model
before it can be trusted. The observed within-autoimmune similarity is compared against 2,000 random gene sets of
the same sizes, drawn from the same gene universe.

![Permutation null versus observed, and the autoimmune-specific genes](/DBRetina/img/autoimmune_proof.png)

The observed within-autoimmune similarity is 40%. The size-matched null sits at 13.6%, and none of the 2,000 random
configurations approached the observed value, which is **167 standard deviations above the null**
(permutation p < 0.001). The clustering is not an artifact of gene-list size.

## What this shows

Six chained DBRetina commands take a set of disease-gene lists to a result that required genome-wide studies to
establish: autoimmune diseases share a genetic core, the core consists of immune and HLA genes, and the
relationship survives a strict null by a wide margin. The same steps apply to any collection of gene sets. The
[reproducible directory](https://github.com/dib-lab/DBRetina) runs the entire workflow with `bash run.sh`.

## References

- Cotsapas C, Voight BF, Rossin E, et al. Pervasive sharing of genetic effects in autoimmune disease.
  *PLoS Genetics* 2011;7(8):e1002254. [doi:10.1371/journal.pgen.1002254](https://doi.org/10.1371/journal.pgen.1002254)
- Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET: a comprehensive platform integrating information on
  human disease-associated genes and variants. *Nucleic Acids Research* 2017;45(D1):D833–D839.
  [doi:10.1093/nar/gkw943](https://doi.org/10.1093/nar/gkw943)
