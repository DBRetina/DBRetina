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


<span class="cmd"> -i, --index-prefix TEXT   index file prefix  [required] </span>

This is the user-defined prefix that was used in the indexing step.

<span class="cmd"> --modularity FLOAT RANGE  containment cutoff for modularity calculation  [default: 80; 0<=x<=100] </span>

This parameter specifies the containment cutoff for modularity calculation. The default value is 80%.

<span class="cmd"> --dedup FLOAT RANGE       deduplication similarity cutoff  [default: 100; 0<=x<=100] </span>

This parameter specifies the Ochiai similarity threshold for deduplication. The default value is 100% which means no items loss and deduplicate only identical groups.

<span class="cmd"> --community FLOAT RANGE   community detection similarity cutoff  [default: 30; 0<=x<=100] </span>

This parameter specifies the Ochiai similarity threshold for community detection. The default value is 30%.

<span class="cmd"> --stop-cov FLOAT RANGE    stop when items covered by %  [0<=x<=100] </span>

This parameter specifies the percentage of items to be covered before stopping the set cover algorithm. The default value is 100% which means the algorithm will stop only when all items are covered.

<span class="cmd"> -o, --output TEXT         output file prefix  [required] </span>

This is the user-defined prefix for the output files.

<hr class="fancy-hr">

## Output files format

<span class="cmd"> {output_prefix}_item_to_GPI_CSI.tsv </span>

A TSV file with the following columns:


<table>
  <tbody>
    <tr>
      <td><strong>Item name</strong></td>
      <td>The item name (i.e. gene, protein)</td>
    </tr>
    <tr>
      <td><strong>GPI</strong></td>
      <td>Group Pleitropy Index</td>
    </tr>
    <tr>
      <td><strong>CSI</strong></td>
      <td>Cluster Specificity Index</td>
    </tr>
  </tbody>
</table>


<span class="cmd"> {output_prefix}_groups_metadata.tsv </span>

A TSV file with the following columns:

<table>
  <tbody>
    <tr>
      <th>Column</th>
      <th>Description</th>
    </tr>
    <tr>
      <td>group</td>
      <td>The group name (node name)</td>
    </tr>
    <tr>
      <td>no_of_items</td>
      <td>Number of items in the group</td>
    </tr>
    <tr>
      <td>average_gpi</td>
      <td>Group's average Pleitropy Index</td>
    </tr>
    <tr>
      <td>average_CSI</td>
      <td>Group's average Cluster Specificity Index</td>
    </tr>
    <tr>
      <td>fragmentation</td>
      <td>Number of outbound edges from the node (-ve value)</td>
    </tr>
    <tr>
      <td>heterogeneity</td>
      <td>Number of inbound edges to the node</td>
    </tr>
    <tr>
      <td>modularity</td>
      <td>absolute(fragmentation + heterogeneity)</td>
    </tr>
    <tr>
      <td>status</td>
      <td>The status of the group.<br>
          `set-cov`: removed by the set-cov.<br> 
          `dedup`: removed by deduplication.<br>
          `remained`: remained group.
      </td>
    </tr>
  </tbody>
</table>

<span class="cmd"> {output_prefix}_remaining_groups_metadata.tsv </span>

Same as {output_prefix}_groups_metadata.tsv but only for the remaining groups after set cover and without the `status` column.

<span class="cmd"> {output_prefix}_removed_groups_metadata.tsv </span>

Same as {output_prefix}_groups_metadata.tsv but only for the removed groups after set cover and without the `status` column.

<span class="cmd"> {output_prefix}_associations.tsv </span>

A new association files as described in [`DBretina index`](dbretina_index.md) section. This file contains the final remaining groups after set cover and deduplication.

<span class="cmd"> {output_prefix}_new.gmt </span>

A new GMT file as described in [`DBretina index`](dbretina_index.md) section. This file contains the final remaining groups after set cover and deduplication.

<span class="cmd"> {output_prefix}_original.gmt </span>

The original GMT file as described in [`DBretina index`](dbretina_index.md) section. This file contains the original groups before set cover and deduplication.

