# Modularity

The `modularity` command calculates the fragmentation, heterogeneity, and modularity of the groups.


```
Usage: DBRetina modularity [OPTIONS]

  Compute the modularity of gene sets

Options:
  -i, --index-prefix TEXT   Index file prefix  [required]
  -p, --pairwise PATH       pairwise TSV file  [required]
  -c, --cutoff FLOAT RANGE  containment cutoff  [0<=x<=100; required]
  -o, --output TEXT         output file prefix  [required]
  --help                    Show this message and exit.
```

## Command arguments

<span class="cmd"> -i, --index-prefix TEXT   Index file prefix  [required] </span>

This is the user-defined prefix that was used in the indexing step as an output prefix.

<span class="cmd">  -p, --pairwise PATH  the pairwise TSV file  [required] </span>

The original or a filtered pairwise TSV file.

<span class="cmd">-c, --cutoff FLOAT RANGE  containment cutoff  [0<=x<=100; required] </span>

The `-c --cutoff` argument uses the **Containment** metric.


<span class="cmd"> -o, --output TEXT    output prefix  [required] </span>

The output files prefix.


<hr class="fancy-hr">


## Output files format

<span class="cmd"> {output_prefix}_modularity.txt </span>

A TSV file that provides information about the fragmentation, heterogeneity, and modularity of the groups. The TSV columns are defined as follows:

<table>
  <tbody>
    <tr>
      <td><strong>group_name</strong></td>
      <td>The group name (node name)</td>
    </tr>
    <tr>
      <td><strong>fragmentation</strong></td>
      <td>Number of outbound edges from the node (-ve value)</td>
    </tr>
    <tr>
      <td><strong>heterogeneity</strong></td>
      <td>Number of inbound edges to the node</td>
    </tr>
    <tr>
      <td><strong>modularity</strong></td>
      <td>absolute(modularity + heterogeneity)</td>
    </tr>
  </tbody>
</table>

The lower modularity value indicates a more modular group. Which means that the group is more likely to be a representative group.