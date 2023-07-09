# Apply set-cover and deduplication algorithms

## Algorithm Overview

DBRetina presents a robust set-cover algorithm that functions within a complex system of parameters. This documentation provides a comprehensive overview of the underlying mechanisms at work, guiding you through the processes.

The set cover algorithm executes several complex computations to achieve its end goals. Here's an ordered breakdown of the process:

1. **Community Detection:** Using the value specified by `--community` (default: 30), the algorithm calculates Ochiai similarity and applies a community detection method on the primary pairwise similarities file, producing clusters of groups. This step serves to calculate the Cluster Specificity Index (CSI) for each group.

2. **Group Pleitropy Index (GPI) Calculation:** Next, the algorithm computes the GPI for each group.

3. **Modularity Calculation:** Using the containment cutoff specified by `--modularity` (default: 80), the algorithm processes the main pairwise similarities file to compute three key metrics: Fragmentation, Heterogeneity, and Modularity. 
   - *Fragmentation* is represented as the negative of the number of outbound edges.
   - *Heterogeneity* is the positive of the number of inbound edges.
   - *Modularity* is the absolute value of the sum of Heterogeneity and Fragmentation. The group with the lowest modularity value has the highest modularity.

4. **Deduplication:** The `--dedup` parameter specifies the Ochiai similarity threshold (default: 100%) for deduplication, creating an undirected graph. The algorithm applies a weakly connected components method, selecting the group with the highest number of edges from each connected component. In case of a tie, the largest group is selected. A default value of 100% is utilized to avoid item loss. Note that the `--dedup` argumnt is similar to the one in [`DBRetina dedup`](dbretina_dedup.md){:target="_blank"}.

5. **Set Cover Application:** Finally, the set cover algorithm sorts the groups by their modularity, CSI, and size, subsequently subtracting shared items from the universal set. This process continues until the `--stop-cov` percentage of items are covered. By default, the algorithm stops only when 100% of items are covered, thereby ensuring maximum coverage.

DBRetina's set cover algorithm offers granular control over the balance between complexity and precision in handling large data sets. By adjusting the various parameters, users can tailor the algorithm's performance to meet specific project requirements.


## Usage


```
Usage: DBRetina setcov [OPTIONS]

  Apply set-cover algorithm.

Options:
  -i, --index-prefix TEXT   index file prefix  [required]
  --modularity FLOAT RANGE  containment cutoff for modularity calculation
                            [default: 80; 0<=x<=100]
  --dedup FLOAT RANGE       deduplication similarity cutoff  [default: 100;
                            0<=x<=100]
  --community FLOAT RANGE   community detection similarity cutoff  [default:
                            30; 0<=x<=100]
  --stop-cov FLOAT RANGE    stop when items covered by %  [0<=x<=100]
  -o, --output TEXT         output file prefix  [required]
  --help                    Show this message and exit.

  Read more at
  https://dbretina.github.io/DBRetina/usage/dbretina_setcov
```

## Command arguments


<span class="cmd"> -i, --index-prefix TEXT  Index file prefix  [required] </span>

This is the user-defined prefix that was used in the indexing step.

<span class="cmd"> -g, --groups-file PATH    single-column supergroups file </span>

This will filter out all pairwise similarities that are between supergroups that are not in the provided groups file. The groups file is a single-column file that contains the names of the supergroups to be included in the filtering.

<span class="cmd"> --cluster-ids TEXT        comma-separated list of cluster IDs </span>

The cluster IDs selected from the clusters file. This argument is only used if the clusters file is not provided.

<span class="cmd"> -o, --output TEXT        output file prefix  [required] </span>

The output prefix that should be unique for this query.

---

## Output files format


<!-- TODO: Implement later -->
<!-- <span class="cmd"> {output_prefix}_features_count_per_group.tsv </span>

A TSV file containing two columns, the first column is the supergroup name and the second column is the number of features that are contained in that supergroup, and the third column is a PIPE-separated list of supergroups from the user query that are associated with the feature in the first column. -->

<span class="cmd"> {output_prefix}_feature_to_groups.tsv </span>

A three-columns TSV file containing the feature name, the number of supergroups associated with that feature, and a PIPE-separated list of supergroups associated with that feature.

<!-- TODO: Fix later -->
<!-- <span class="cmd"> {output_prefix}_features_count_per_group_histogram.png </span>

A histogram plot showing the distribution of the number of features per supergroup. -->