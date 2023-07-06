# 4. Cluster

Graph-based clustering of the pairwise TSV file based on the connected components algorithm. The clustering is based on the distance metric and the cutoff value. The output is a DBRetina clusters file.

```
Usage: DBRetina cluster [OPTIONS]

  Graph-based clustering of the pairwise TSV file.

Options:
  -p, --pairwise PATH       filtered pairwise TSV file
  -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont',
                            'ochiai', 'jaccard']  [required]
  -c, --cutoff FLOAT RANGE  cluster the supergroups with (distance > cutoff)
                            [default: 0.0; 0<=x<=100]
  -o, --output-prefix TEXT  output file prefix  [required]
  --help                    Show this message and exit.
```


## 4.1 Command arguments

<span class="cmd"> -c, --cutoff FLOAT RANGE  cluster the supergroups with (distance > cutoff) [default: 0.0; 0<=x<=100] </span>

The cutoff value for clustering the supergroups. The default value is 0.0, which means that all supergroups will be clustered together. The cutoff value can be between 0% and 100%.

<span class="cmd"> -p, --pairwise PATH       filtered pairwise TSV file  [required] </span>

The original or a filtered pairwise TSV file.

<span class="cmd"> -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']  [required] </span>

The distance metric to apply the cutoff on.

---

## 4.2 Output files format

<span class="cmd"> {output_prefix}_clusters.tsv </span>

The DBRetina clusters TSV file. First column is the cluster ID, second column is the cluster size, and the third column is PIPE separated cluster members.


<span class="cmd"> {output_prefix}_clusters_histogram.png </span>

A histogram provides a visual representation of the distribution of cluster sizes. Each bar corresponds to a size range, with the height of the bar indicating the number of clusters falling within that range. This allows for a quick understanding of how cluster sizes are distributed, identifying common sizes and outliers.

<span class="cmd"> {output_prefix}_clusters_bubbles.png </span>

Bubble plot uses a grid layout to represent distinct clusters. The bubble size and color gradient both denote the magnitude of each cluster. The bubble plot is useful for visualizing the distribution of cluster sizes and the relative sizes of each cluster.
