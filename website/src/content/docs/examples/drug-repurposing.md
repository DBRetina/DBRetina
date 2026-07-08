---
title: "Recovering drug mechanism from shared gene targets"
description: Grouping cancer drugs by their gene targets, and reading off repurposing hints.
---

Most cancer drugs hit more than one gene. If a drug's mechanism is written in its targets, then drugs with the
same mechanism should share targets, and drugs that share targets across mechanisms are repurposing candidates.
The question is how much of that signal the targets actually carry.

The panel is 24 cancer drugs and their target genes, drawn from DSigDB (Yoo et al., 2015), spanning several known
classes:
BCR-ABL inhibitors (imatinib, dasatinib, nilotinib, ponatinib), EGFR inhibitors (gefitinib, erlotinib,
lapatinib), VEGFR and multi-kinase inhibitors (sorafenib, sunitinib, cabozantinib), BRAF inhibitors (vemurafenib,
dabrafenib), and a few drugs with unrelated mechanisms (the HDAC inhibitor vorinostat, the antimetabolite
methotrexate, the anti-estrogen tamoxifen).

```bash
DBRetina index -g drugs.gmt -o drugs
DBRetina pairwise -i drugs -m ochiai -c 0
```

![Clustered heatmap of 24 cancer drugs by shared gene targets, coloured by mechanism](/DBRetina/img/drug_clustermap.png)

The colour bar encodes the mechanism labels, which are not used in the clustering; the groupings recover them.
The EGFR inhibitors form the tightest group, with erlotinib and gefitinib sharing 75% of their targets. The VEGFR
and multi-kinase inhibitors collapse into one broad block together with several BCR-ABL inhibitors, because those
drugs are deliberately promiscuous and hit overlapping kinase families. The three drugs with unrelated
mechanisms, vorinostat, methotrexate, and tamoxifen, sit off on their own with almost no similarity to the rest.

## A repurposing hint

`neighbors` ranks the drugs that share sorafenib's targets:

```bash
DBRetina neighbors -d drugs_DBRetina_pairwise "sorafenib" -m ochiai -c 0
```

```text
neighbor    ochiai  jaccard  shared_features
imatinib    56.1    38.8     97
erlotinib   51.9    33.7     65
gefitinib   51.4    33.2     64
lapatinib   40.5    21.5     38
```

Sorafenib (approved for kidney and liver cancer) shares 97 target genes with imatinib (approved for CML and GIST),
and around 65 with the EGFR inhibitors. Shared targets across indications like this are the signal that
target-based repurposing screens look for.

Target overlap carries real mechanistic signal. The known drug classes emerge from targets alone, and the
cross-class overlaps worth a second look surface in the same comparison. `neighbors` turns that into a shortlist
for any drug.

## References

- Yoo M, Shin J, Kim J, et al. DSigDB: drug signatures database for gene set analysis. *Bioinformatics*
  2015;31(18):3069–3071. [doi:10.1093/bioinformatics/btv313](https://doi.org/10.1093/bioinformatics/btv313)
