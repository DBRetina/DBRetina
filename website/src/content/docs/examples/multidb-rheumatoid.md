---
title: "From disease to mechanism to drug: rheumatoid arthritis across five databases"
description: Combining a disease, a drug, and three pathway databases to find the mechanisms and drugs behind rheumatoid arthritis.
---

A single database only answers questions about its own contents. A disease-gene database tells you which diseases
resemble each other; a drug-target database tells you which drugs resemble each other; neither can connect a
disease to its drugs or its pathways. Because DBRetina indexes any gene sets together, several databases can be loaded into one network and queried
at once.

Five databases are combined into a single index: DisGeNET for diseases (Piñero et al., 2017), DSigDB for drug
targets (Yoo et al., 2015), and KEGG (Kanehisa & Goto, 2000), Reactome (Gillespie et al., 2022), and WikiPathways
(Martens et al., 2021) for pathways. That is about 2,600 gene sets of three different kinds, queried for what
connects to rheumatoid arthritis.

```bash
DBRetina index -g diseases.gmt -g drugs.gmt -g kegg.gmt -g reactome.gmt -g wikipathways.gmt -o multidb
DBRetina pairwise -i multidb -m containment -c 0
DBRetina neighbors -d multidb_DBRetina_pairwise "dis_rheumatoid_arthritis" -m containment -c 0
```

Containment is used here rather than ochiai, because the comparison is directional: how much of a pathway or a
drug's genes fall inside rheumatoid arthritis's much larger gene set. The neighbours split cleanly into the two
other kinds of database:

![Rheumatoid arthritis's top pathways and drugs across five databases](/DBRetina/img/multidb_ra.png)

On the mechanism side, the pathways most contained in rheumatoid arthritis are immune and inflammatory:
interleukin-21 signalling and STAT5 activation from Reactome, the CLEC7A inflammasome, macrophage markers, and
T-cell chemokine polarisation from WikiPathways, and, as a sanity check, a WikiPathways set literally named
"genes associated with the development of rheumatoid arthritis." On the treatment side, the drugs are the ones
rheumatologists prescribe: prednisone, hydroxychloroquine, sulfasalazine, azathioprine, leflunomide, tofacitinib,
and methotrexate, the standard disease-modifying antirheumatic drugs.

The mechanism and the treatment of a disease are assembled from three different kinds of database, using only the
genes they have in common. The same query pointed at another disease returns that disease's pathways and drugs.

## References

- Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET: a comprehensive platform integrating information on
  human disease-associated genes and variants. *Nucleic Acids Research* 2017;45(D1):D833–D839.
  [doi:10.1093/nar/gkw943](https://doi.org/10.1093/nar/gkw943)
- Yoo M, Shin J, Kim J, et al. DSigDB: drug signatures database for gene set analysis. *Bioinformatics*
  2015;31(18):3069–3071. [doi:10.1093/bioinformatics/btv313](https://doi.org/10.1093/bioinformatics/btv313)
- Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 2000;28(1):27–30.
  [doi:10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway knowledgebase 2022. *Nucleic Acids Research*
  2022;50(D1):D687–D692. [doi:10.1093/nar/gkab1028](https://doi.org/10.1093/nar/gkab1028)
- Martens M, Ammar A, Riutta A, et al. WikiPathways: connecting communities. *Nucleic Acids Research*
  2021;49(D1):D613–D621. [doi:10.1093/nar/gkaa1024](https://doi.org/10.1093/nar/gkaa1024)
