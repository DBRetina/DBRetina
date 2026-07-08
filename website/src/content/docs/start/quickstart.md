---
title: Quickstart
description: Build an index, compare every pair, and cluster the results, in about five minutes.
---

This walkthrough takes you from a raw gene-set file to a clustered similarity network. It uses a tiny example so
you can run every command yourself and see the same numbers.

Grab the example file (six diseases and their genes) and follow along:
[example.gmt](/DBRetina/example.gmt).

```text title="example.gmt"
breast_cancer      na  BRCA1 BRCA2 TP53 PTEN ATM CHEK2 PALB2
ovarian_cancer     na  BRCA1 BRCA2 TP53 RAD51C BRIP1 PALB2
lung_cancer        na  EGFR KRAS TP53 ALK MET ROS1
melanoma           na  BRAF NRAS TP53 CDKN2A PTEN KIT
alzheimer_disease  na  APP PSEN1 PSEN2 APOE MAPT TREM2
parkinson_disease  na  SNCA LRRK2 PARK7 PINK1 PRKN GBA
```

## 1. Build an index

The index is DBRetina's internal representation of your groups and their features. Build it once; every other
command reads it.

```bash
DBRetina index -g example.gmt -o example
```

This writes three files: `example.dbri` (the index) plus `example_raw.json` and `example_hashes.json` (your
groups in plain text and hashed form).

## 2. Compare every pair

```bash
DBRetina pairwise -i example -m ochiai -c 0
```

`pairwise` compares all groups against each other and writes the results to a directory,
`example_DBRetina_pairwise/`. Every pair gets a row with the number of shared features and six similarity
metrics. Here `-m ochiai -c 0` keeps all pairs (cutoff of 0); raise the cutoff to keep only the strong ones.

## 3. Look at what's similar

The quickest way to read the result is to ask for a group's nearest neighbours:

```bash
DBRetina neighbors -d example_DBRetina_pairwise "breast_cancer" -m ochiai -c 0
```

```text
[INFO] Showing top 3 neighbors of "breast_cancer" (ochiai >= 0.0)
neighbor         ochiai  jaccard  shared_features
ovarian_cancer   61.7    44.4     4
melanoma         30.9    18.2     2
lung_cancer      15.4    8.3      1
```

Breast and ovarian cancer share four genes (BRCA1, BRCA2, TP53, PALB2), so they come out most similar, which is
exactly what you'd expect.

## 4. Cluster the network

Group the diseases by their similarity. Community detection at an ochiai cutoff of 20 finds the connected family:

```bash
DBRetina cluster -p example_DBRetina_pairwise -m ochiai -c 20 --community -o example_clusters
```

```text
[INFO] Total number of clustered supergroups: 3
[INFO] number of clusters: 1
```

```text title="example_clusters_clusters.tsv"
cluster_id  cluster_size  cluster_members
1           3             breast_cancer|ovarian_cancer|melanoma
```

The three cancers cluster together (they share TP53 and PTEN); the two neurodegenerative diseases share nothing
with them and drop out at this cutoff.

## 5. Draw the tree

To see the whole set as a dendrogram and heatmap:

```bash
DBRetina export -p example_DBRetina_pairwise -m ochiai --newick -o example_tree
```

This writes `example_tree_dendrogram.png`, `example_tree_heatmap.png`, a distance matrix, and a Newick tree you
can open in any tree viewer.

## Where to go next

- Read [Core concepts](/DBRetina/start/concepts/) to understand the metrics and the files.
- Work through a real, multi-database analysis in [Examples](/DBRetina/examples/).
- Look up any command's exact options under [Commands](/DBRetina/commands/).
