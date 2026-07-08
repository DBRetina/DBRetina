---
title: "Comorbidity, not cause: autism's genetic neighbours"
description: Recovering autism's coexisting neuropsychiatric conditions from gene-set overlap, and separating polygenic comorbidity from single-gene cause.
---

Autism co-occurs with a broad range of neuropsychiatric conditions. Epilepsy, ADHD, and intellectual disability
are common in childhood, and the risk of anxiety, depression, bipolar disorder, and schizophrenia is elevated into
adulthood. This comorbidity has a genetic component: large psychiatric genetics studies have shown that these
disorders share genetic risk (Cross-Disorder Group of the Psychiatric Genomics Consortium, 2013). The analysis
that follows tests whether that structure is recoverable from disease-gene lists alone, and whether the same
gene-set overlap distinguishes the polygenic comorbidities from the single-gene syndromes (Rett, Fragile X,
tuberous sclerosis) that cause autism.

The panel comprises 18 disorders drawn from DisGeNET (Piñero et al., 2017): autism itself; the neurodevelopmental
conditions it most frequently co-occurs with (epilepsy, ADHD, intellectual disability, developmental delay); the
psychiatric disorders (schizophrenia, bipolar, depression, anxiety, OCD, Tourette); the monogenic syndromes (Rett,
Fragile X, tuberous sclerosis, Down); and three controls drawn from outside neurology (rheumatoid arthritis, breast
cancer, coronary artery disease).

## Building the network

`index` hashes the gene sets into DBRetina's internal representation, `pairwise` computes all-versus-all ochiai
similarity with the cutoff at zero so that no edges are dropped, and `export` renders the resulting matrix as a
dendrogram.

```bash
DBRetina index -g neuro.gmt -o nz
DBRetina pairwise -i nz -m ochiai -c 0
DBRetina export -p nz_DBRetina_pairwise -m ochiai --newick -o nz_tree -l names
```

![Circular dendrogram of the 18 disorders](/DBRetina/img/neuro_dendrogram.png)

The dendrogram is computed from gene-set overlap alone, with no group labels supplied, and it recovers the clinical
groupings. Autism falls next to epilepsy and ADHD, the neurodevelopmental triad with the highest clinical
co-occurrence. The mood and psychotic disorders (schizophrenia, bipolar, depression, anxiety) form a second clade;
the single-gene syndromes (Rett, Fragile X, tuberous sclerosis) a third; and intellectual disability pairs with
developmental delay. The three controls group together, separate from the neurological disorders.

A dendrogram imposes a hierarchy and obscures the direct pairwise links. The `graph` command represents the same
similarity matrix as a node-and-edge network. Retaining edges above 22% similarity and colouring each node by its
group:

```bash
DBRetina graph -i nz -p nz_DBRetina_pairwise -m ochiai -c 22 -o nz_graph
```

![Force-directed network of autism and its coexisting disorders](/DBRetina/img/autism_graph.png)

Autism occupies a central position, connected to the neurodevelopmental disorders on one side and the psychiatric
disorders on the other. The controls separate toward the top right. The single-gene syndromes attach at the
periphery through low-weight edges.

## Autism's nearest neighbours

`module` extracts the disorders most similar to a query group and returns them as a set. For autism:

```bash
DBRetina module -d nz_DBRetina_pairwise dis_asd -m ochiai -c 10 --min-shared 20 -n 6
```

```text
schizophrenia
epilepsy
adhd
anxiety
intellectual_disability
depression
```

All six are documented autism comorbidities.

## Comorbidity, not cause

The neighbour set is the top of a full ranking. Plotting autism's similarity to every disorder in the panel,
coloured by group, resolves three bands.

![Autism's similarity to every disorder in the panel, coloured by group](/DBRetina/img/autism_neighbours.png)

Autism's seven strongest neighbours are all polygenic psychiatric or neurodevelopmental comorbidities:
schizophrenia, epilepsy, ADHD, anxiety, intellectual disability, depression, and bipolar. Each ranks above every
non-neuro control (Mann-Whitney p = 0.008).

This separation is specific to autism rather than a property of the neurological group as a whole. Autism's
neighbours clear the controls, but the neurological disorders taken together do not, because the monogenic
syndromes fall below the controls. Rett, Fragile X, and tuberous sclerosis average 18.4% similarity to autism,
among the lowest values in the panel, controls included, despite each of the three causing autism.

The inversion follows from gene-set architecture. Each syndrome is driven by a single gene (MECP2 in Rett, FMR1 in
Fragile X, TSC2 in tuberous sclerosis), so its DisGeNET gene set is small and narrowly defined. Idiopathic autism
is polygenic, distributed across hundreds of genes, and its overlap aligns with the other polygenic disorders that
share this architecture rather than with the narrow single-gene definitions. The distinction between comorbidity
and cause is encoded in the structure of the gene sets themselves.

## The shared gene core

`geneinfo` reports, for a specified set of disorders, the genes common to them. Applied to autism and its core
comorbidities:

```bash
DBRetina geneinfo -i nz -g core.txt -o core_genes   # core.txt = autism, epilepsy, ADHD, ID, schizophrenia
```

The shared genes are synaptic and neurodevelopmental, the genes recurrently implicated by psychiatric exome
studies: SHANK3, GRIN2B, GRIN2A, CNTNAP2, SYNGAP1, NRXN1, GABRB3, RELN, MEF2C, and the SCN-family channels
(Satterstrom et al., 2020). The three single-gene culprits, MECP2, FMR1, and TSC2, are present in this core as
well. The genes belong to the broad polygenic disorders and appear in their overlap, whereas the syndromes named
for them are defined narrowly enough to fall outside the cluster.

## What this shows

From disease-gene lists, the analysis recovers autism's comorbidity structure, identifies a shared genetic core
consistent with the architecture established by genome-wide psychiatric genomics, and separates the polygenic
comorbidities from the monogenic syndromes that cause autism. The full workflow, including commands and inputs, is
available in the [reproducible directory](https://github.com/dib-lab/DBRetina).

## References

- Cross-Disorder Group of the Psychiatric Genomics Consortium. Genetic relationship between five psychiatric
  disorders estimated from genome-wide SNPs. *Nature Genetics* 2013;45(9):984–994.
  [doi:10.1038/ng.2711](https://doi.org/10.1038/ng.2711)
- Satterstrom FK, Kosmicki JA, Wang J, et al. Large-scale exome sequencing study implicates both developmental
  and functional changes in the neurobiology of autism. *Cell* 2020;180(3):568–584.
  [doi:10.1016/j.cell.2019.12.036](https://doi.org/10.1016/j.cell.2019.12.036)
- Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET. *Nucleic Acids Research* 2017;45(D1):D833–D839.
  [doi:10.1093/nar/gkw943](https://doi.org/10.1093/nar/gkw943)
