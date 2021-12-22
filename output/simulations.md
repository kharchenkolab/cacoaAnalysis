# Simulations for evaluation of expression change metrics

In these simulations we used Muscat **REF** for simulation of artificial datasets with parameters estimated based on the Autism data. We varied different parameters, such as number of cells and samples and compared several metrics of gene expression change:

- correlation distances between expression shifts ("raw distances")
- cluster-based expression shifts
- cluster-based expression shifts estimated on top-500 differentially expressed (DE) genes with subsequent PCA
- cluster-free expression shifts
- number of significant DE genes

In most cases, we generated data based on one cell type ("IN-PV") with different parameters, 30 times for each combination of parameters.

## Dependency on the number of cells and samples

### Number of samples per cell type

We start by showing the raw expression distances (y-axis) depending on the number of samples (x-axis). Hereafter, each dot represent median distance within one simulations. We can see that raw distances do not depend on the number of samples:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-3-1.png" width="480" />

The normalized distance does not depend either, however the statisticsl power of the test depends a lot:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-4-1.png" width="864" />

Estimating expression shifts over top DE introduce slight dependency. Probably, because of the quality of the selected DE genes.

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-5-1.png" width="864" />

Cluster-free estimates also do not have the dependency:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-8-1.png" width="768" />

The number of significant DE genes strongly depends on the number of samples:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-9-1.png" width="480" />

It happens because the number of detected genes depends on the power of the test, which depends on the number of samples.

### Number of cells per cell type

The raw expression distances depend on the number of samples:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-11-1.png" width="480" />

However the normalization fixes it:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-12-1.png" width="864" />

And here, using top-DE genes does not introduce new dependencies:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-13-1.png" width="864" />

The same with cluster-free:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-15-1.png" width="768" />

Jsut as in the previous case, the number of significant DE genes does depend on the number of cells:

<img src="../docs/figure/simulation_ns_nc.Rmd/unnamed-chunk-16-1.png" width="480" />

## Simulations for different cell types

### Dependency on the number of genes and other cell type covariates

Muscat doesn’t allow varying number of genes per cell type properly, so we simulated data from different cell types and visualized dependencies of the resulting distances on the parameters of these cell types. It does not allow to analyze possible sources of variation independently, but still provides some understanding of the dependencies.

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-4-1.png" width="576" />

Here, we additionally varied fraction of simulated DE genes (color) to compare how much of signal comes from it compared to the signal from the covariates.

It can be seen that there is a lot of variation for the same DE fraction, even if it is set to 0.0. To explain the variation we may plot different cell type specific covariates.

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-5-1.png" width="672" />

This plot shows the number of expressed genes in real data (bottom), mean number of molecules per gene in real data (middle) and median distance between samples within the same condition for DE fraction set to 0.00 (top).

First, it can be seen that the distances within the same cell type are generally higher for cell types with low coverage. Second, variation between samples within the condition appear to be the main driver of the distances between the conditions.

Normalizing distances as we do to estimate expression shifts we can greatly reduce those biases:

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-6-1.png" width="864" />

And here, focusing on top genes reduces biases even more:

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-7-1.png" width="864" />

Cluster-free estimates are the most sensitive to variation in the data. They still work to prioritize cell types with adequate variation, but are likely to miss variable cell types:

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-9-1.png" width="768" />

All the plots above show that the distances for DE fraction set to 0.0 are slightly below 0.0. It is an artifact of simulation, but not one of the distances. This happens because with zero DE genes, Muscat generates two identical sets of samples. So, when estimating distances **between** conditions we have exactly the same distances as we have when estimating **within** the same conditions with the addition of one value of 0.0 as the distance to the identical sample. So, the distances between conditions are on average **lower** than distances within the same condition.

Finally, below is the same plot for the number of significant DE genes. By design of the simulations, the number of DE genes linearly depends on the total number of genes. So, the right plot shows the fraction of DE genes, which also depends (slightly) on the variation within each cell types. It is up to debate, which of these measures should be used for real world examples, though.

<img src="../docs/figure/simulation_types.Rmd/unnamed-chunk-10-1.png" width="960" />

## Sensitivity to expression change

The two section below test how much different measures change when actual change in the expression happens. We mesure amount of change as log2-fold change (LFC) and the fraction of DE genes.

### Sensitivity to log2-fold change (LFC)

Muscat does not allow simulating data with LFC &lt; 1, so we analyse only values above 1. DE fraction is fixed to 0.05.

Increasing log2-fold change affects the normalized distance a lot:

<img src="../docs/figure/simulation_sensitivity_lfc.Rmd/unnamed-chunk-3-1.png" width="864" />

Using top DE genes does not increase sensitivity here, as it is already high for DE fraction of 0.05 and LFC > 1.

<img src="../docs/figure/simulation_sensitivity_lfc.Rmd/unnamed-chunk-4-1.png" width="864" />

Cluster-free estimates have lower sensitivity, but still reach significance pretty fast:

<img src="../docs/figure/simulation_sensitivity_lfc.Rmd/unnamed-chunk-6-1.png" width="768" />

The number of DE genes also gets more sensitive as the LFC increases:

<img src="../docs/figure/simulation_sensitivity_lfc.Rmd/unnamed-chunk-7-1.png" width="480" />


### Sensitivity to DE fraction

The normalized distances depend linearly on the fraction of DE genes and reach significance pretty fast:

<img src="../docs/figure/simulation_sensitivity_frac.Rmd/unnamed-chunk-3-1.png" width="864" />

Here, we can see that using top DE genes increases sensitivity for low number of DE genes. However, after certain amount of changes, this distance reaches plateau, because adding more DE genes does not contribute to the amount of changes in top.

<img src="../docs/figure/simulation_sensitivity_frac.Rmd/unnamed-chunk-4-1.png" width="864" />

Cluster-free shifts are much less sensitive and more variable than the cluster-based version:

<img src="../docs/figure/simulation_sensitivity_frac.Rmd/unnamed-chunk-6-1.png" width="768" />

The number of DE genes here behaves exactly as expected:

<img src="../docs/figure/simulation_sensitivity_frac.Rmd/unnamed-chunk-7-1.png" width="480" />

## Evaluation of different normalization methods for expression distances

These simulations are very different from the previous cases. Here we did not use Muscat. Insted we artificially generated pairwise distance matrix by sampling from two normal distributions with controlled distance between their means and controlled variance of each distribution. Then, we compared three types of distance estimates:

1.  **Shift**: <span class="math inline">\\(\\textbf{d}\_{between} - 0.5 \\cdot \\left( \\mathrm{MED}(\\textbf{d}\_{case}) + \\mathrm{MED}(\\textbf{d}\_{control}) \\right)\\)</span>
2.  **Total**: <span class="math inline">\\(\\textbf{d}\_{between} - \\mathrm{MED}(\\textbf{d}\_{control})\\)</span>
3.  **Variance**: <span class="math inline">\\(\\textbf{d}\_{case} - \\mathrm{MED}(\\textbf{d}\_{control})\\)</span>

### Equal variance, varying mean

First, we changed mean between two conditions, keeping the variance constant.

<img src="../docs/figure/simulation_distances.Rmd/unnamed-chunk-1-1.png" width="576" />

The plot shows change in expression distances (y-axis) of different types (color) for increasing mean offset between case and control (x-axis) and fixed standard deviation (std=1.0 for the top plot and std=2.0 for the bottom plot). It can be seen that the "Variance" distance does depend on mean offset, while "Shift" and "Total" go up linearly. It is also important to note that when variance within condition increase (bottom plot), the sensitivity of all distances goes down. So, if two cell types have the same mean shift, but different variance, the cell type with lower variance would have higher normalized distance.

### Equal mean, varying variance

Now, we will fix the mean offset to different values and vary variance within Case samples, keeping Control variance fixed (`std=1.0`).

<img src="../docs/figure/simulation_distances.Rmd/unnamed-chunk-2-1.png" width="576" />

This panel shows change in expression distances (y-axis) of different types (color) for increasing standard deviation in Case (x-axis) and fixed Mean offset (from 0 to 2 in the three plots top-down). It can be seen that the "Shift" distance is quite insensitive to the changes in variance, while "Variance" distance behaving exactly the same for each offset value. It also shows that "Total" distance captures both changes.
