# 2. Pairwise

The Pairwise command in DBRetina is designed to perform pairwise comparisons between supergroups based on their shared features. This command takes the index prefix and the number of cores as input parameters.


```
Usage: DBRetina pairwise [OPTIONS]

  Calculate pairwise distances.

Options:
  -i, --index-prefix TEXT   Index file prefix  [required]
  -t, --threads INTEGER     number of cores
  -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont',
                            'ochiai', 'jaccard']  [default: max_cont]
  -c, --cutoff FLOAT RANGE  filter out distances < cutoff  [default: 0.0;
                            0<=x<=100]
  --help                    Show this message and exit.
```

## 2.1 Command arguments

<span style="color:orange;">** -i, --index-prefix TEXT  Index file prefix  [required] **</span>

This is the user-defined prefix that was used in the indexing step.

<span style="color:orange;">** -t, --threads INTEGER    number of cores **</span>

The number of processing cores to be used for parallel computation during the pairwise comparisons.

<span style="color:orange;">** -d, --dist-type TEXT     select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']  [default: max_cont] **</span>

<span style="color:orange;">** -c, --cutoff FLOAT RANGE filter out distances < cutoff  [default: 0.0; 0<=x<=100] **</span>

The `-d` and `-c` input parameters serve the purpose of selecting a particular distance metric and predefined cutoff. This cutoff will eliminate all pairwise comparisons that have a distance value lower than the cutoff.

---

## 2.2 Output files format

<span style="color:orange;">** {perfix}_DBRetina_pairwise.tsv **</span>

A TSV file that provides information about shared features between each pair of supergroups. The TSV columns are defined as follows:


<table>
  <tbody>
    <tr>
      <td><strong>group_1_ID</strong></td>
      <td>ID of the first supergroup in a pair</td>
    </tr>
    <tr>
      <td><strong>group_2_ID</strong></td>
      <td>ID of the second supergroup in a pair</td>
    </tr>
    <tr>
      <td><strong>group_1_name</strong></td>
      <td>name of the first supergroup in a pair</td>
    </tr>
    <tr>
      <td><strong>group_2_name</strong></td>
      <td>name of the second supergroup in a pair</td>
    </tr>
    <tr>
      <td><strong>shared_features</strong></td>
      <td>number of features shared between the two supergroups</td>
    </tr>
    <tr>
      <td><strong>min_containment</strong></td>
      <td>minimum containment between the two supergroups</td>
    </tr>
    <tr>
      <td><strong>avg_containment</strong></td>
      <td>average containment between the two supergroups</td>
    </tr>
    <tr>
      <td><strong>max_containment</strong></td>
      <td>maximum containment between the two supergroups</td>
    </tr>
    <tr>
      <td><strong>ochiai</strong></td>
      <td>Ochiai distance between the two supergroups</td>
    </tr>
    <tr>
      <td><strong>jaccard</strong></td>
      <td>Jaccard distance between the two supergroups</td>
    </tr>
  </tbody>
</table>


# The output PNG file of histogram of pairwise distances

<span style="color:orange;">** {index_prefix}_DBRetina_distance_metrics_plot_log.png **</span>

clustered bar chart illustrates the frequency distribution of five distance metrics - min_cont, avg_cont, max_cont, ochiai, and jaccard - across various distance ranges. The y-axis is displayed on a logarithmic scale to accommodate the wide range of frequencies observed in the data.

<span style="color:orange;">** {index_prefix}_DBRetina_distance_metrics_plot_linear.png **</span>

Same as above, but the y-axis is displayed on a linear scale.