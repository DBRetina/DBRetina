---
title: Metrics
description: The similarity metrics DBRetina computes, and how to read them.
---

For a pair of groups, let `s1` and `s2` be their feature-set sizes and `o` the number of features they share.
DBRetina reports the following, each scaled to a **0–100 percentage**.

## Shared features

```text
shared_features = o
```

The raw count of features in common. Everything else is a normalised version of this.

## Containment

```text
containment = o / min(s1, s2)
```

How much of the *smaller* group sits inside the larger one. Containment is asymmetric in spirit: it reaches 100
when one group is a subset of the other, which makes it the right choice for "is this group part of that one".

## Ochiai

```text
ochiai = o / sqrt(s1 · s2)
```

The geometric mean of the two containments. It balances the two group sizes, so a small group and a large group
only score high if the overlap is substantial relative to both. A good general-purpose symmetric metric.

## Jaccard

```text
jaccard = o / (s1 + s2 − o)
```

Overlap over the union. The strictest of the symmetric metrics: two groups score high only when they are alike in
both membership and size.

## CSI

```text
csi = o² / (s1 · s2)
```

The squared cosine similarity, that is Ochiai squared, which sharpens the gap between strong and weak overlaps.
Despite the abbreviation, this is **not** the graph-theoretic Connection Specificity Index; it is a simple
size-aware overlap.

## Dice

```text
dice = 2o / (s1 + s2)
```

Overlap over the average size. Similar in spirit to Jaccard but more forgiving of size differences.

## Odds ratio

```text
odds_ratio = (o · d) / ((s1 − o) · (s2 − o)),   d = N − (s1 + s2 − o)
```

where `N` is the feature universe, all distinct features across the database. This is the odds ratio of the 2×2
table of features in both groups, in one only, or in neither: how much more the two groups overlap than if their
features were drawn at random from the universe. Above 1 means more overlap than chance. It is undefined (reported
as −1 and dropped) when one group's features are wholly contained in the other, so unlike the ratios above it is a
size-aware significance score, not a bounded 0–100 similarity.

## P-value (optional)

With `pairwise --pvalue`, DBRetina adds a hypergeometric p-value testing whether the observed overlap is larger
than chance. Unlike every metric above, **lower is better**: a small p-value means a surprising overlap. A
p-value cutoff therefore keeps the pairs at or below the threshold, and `--fdr` adds Benjamini–Hochberg
correction across all pairs.

## Choosing one

- **Containment**: subset relationships, "is A inside B".
- **Ochiai**: balanced, size-aware similarity; a sensible default.
- **Jaccard / Dice**: strict symmetric similarity when size should matter.
- **Odds ratio / p-value**: statistical significance rather than raw similarity.
