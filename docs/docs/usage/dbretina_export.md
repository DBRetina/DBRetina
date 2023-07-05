# 4. Export

The Export command converts the pairwise TSV file to a dissimilarity matrix and (optionally) a newick-format file.

```
Usage: DBRetina export [OPTIONS]

  Export to dissimilarity matrix and newick format.

  Export a pairwise TSV file to a dissimilarity matrix and (optionally) a
  newick-format file.

Options:
  -p, --pairwise PATH   filtered pairwise TSV file  [required]
  -d, --dist-type TEXT  select from ['min_cont', 'avg_cont', 'max_cont',
                        'ochiai', 'jaccard']  [default: max_cont]
  --newick              Convert the dissimilarity matrix to newick tree format
  -l, --labels TEXT     select from ['ids', 'names']  [default: ids]
  -o TEXT               output prefix  [required]
  --help                Show this message and exit.
```

## 4.1 Command arguments


<span style="color:orange;">** -p, --pairwise PATH       the pairwise TSV file  [required] **</span>

The original or a filtered pairwise TSV file.

<span style="color:orange;">** -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']  [default: max_cont] **</span>

The distance metric to be utilized in the similarity matrix.

<span style="color:orange;">** --newick              Convert the similarity matrix to newick tree format **</span>

This will convert the similarity matrix to a newick tree format.

<span style="color:orange;">** -l, --labels TEXT         select from ['ids', 'names']  [default: ids] **</span>

The labels to be used in all outputs (plots, newick, distance matrix) matrix. The default is the IDs which are the IDs found in the {index_prefix}.namesMap file.

<span style="color:orange;">** -o TEXT               custom output file name prefix **</span>

The output files prefix.


---


## 4.2 Output files format

<span style="color:orange;">** {output_prefix}_distmat.tsv **</span>

A labels tab-separated similarity matrix file.

<span style="color:orange;">** {output_prefix}_distmat.pkl **</span>

A labels pickle similarity matrix file in Python pickle format.

<span style="color:orange;">** {output_prefix}.newick **</span>

A newick tree file that can be used in tree visualization tools.

<span style="color:orange;">** {output_prefix}_clustermap.png **</span>

A clustermap plot of the similarity matrix.