Alzheimer’s disease, entorhinal cortex: the cluster-free maps on a
Harmony-aligned graph
================
Cacoa developers
2026-11-15

- <a href="#load-data" id="toc-load-data">Load data</a>
- <a href="#cluster-free-density-shifts"
  id="toc-cluster-free-density-shifts">Cluster-free density shifts</a>
- <a href="#cluster-free-expression-shifts"
  id="toc-cluster-free-expression-shifts">Cluster-free expression
  shifts</a>
- <a href="#session-info" id="toc-session-info">Session info</a>

This notebook is a companion of `alz_cacoa.Rmd`: the same cells, counts
and annotation on a kNN graph of a Harmony embedding
(`prepare_alz_harmony.R`: 2,000 variable genes, 50 PCs, Harmony by
sample, 30 nearest neighbors on 30 Harmony dimensions, UMAP on the same)
instead of the Conos alignment graph. Only the graph-dependent parts are
rerun, the cluster-free density and expression-shift maps; the
cell-type-level analyses use the counts and the annotation alone and are
those of the main notebook. The PF companion
(`../pf/pf_cacoa_harmony.Rmd`, `../pf/pf_density_diagnosis.md`) shows
why the graph matters there: on the Conos alignment of that cohort the
density kernel leaks out of the populations present in only part of the
samples. Here the question is whether the maps of a cohort where the
Conos graph worked change with the alignment.

``` r
library(cacoa)
library(ggplot2)
library(dplyr)
```

## Load data

``` r
seu <- readRDS(params$data.file)
cell.meta <- read.delim(params$meta.file)
sample.meta <- cell.meta %>% filter(patient %in% unique(seu$sample)) %>% group_by(patient) %>%
  summarise(condition = ifelse(unique(batchCond) == "AD", "AD", "Control"), sex = unique(sex), batch = unique(batch),
            n.cells = n(), median.umi = median(nUMI), median.genes = median(nGene), median.pct.mt = median(pctMT)) %>%
  mutate(condition = factor(condition, levels = c("Control", "AD"))) %>% as.data.frame()
rownames(sample.meta) <- sample.meta$patient
sample.meta$patient <- NULL
sample.meta <- sample.meta[sort(unique(seu$sample)), ]
table(sample.meta$condition)
```

    ## 
    ## Control      AD 
    ##       6       6

``` r
cao <- Cacoa$new(seu, sample.metadata = sample.meta, sample.per.cell = factor(seu$sample), cell.groups = factor(seu$celltype),
                 embedding = seu@reductions$umap@cell.embeddings, graph.name = "RNA_nn", data.layer = "data", n.cores = params$n.cores)
```

    ## Sample metadata: 12 samples, 7 columns
    ##   usable (6): condition [AD 6, Control 6], sex [Female 4, Male 8], n.cells [num, 305-1964], median.umi [num, 542-1234], median.genes [num, 434-824], median.pct.mt [num, 0.0562-0.923]
    ##   excluded from screening (1): batch (high-cardinality)

    ## Next: cao$checkDesign() and cao$screenCovariates(), or cao$setModel(test = "<variable>").

``` r
cao$setOptions(n.permutations = 999, seed = 1)
cao$plot.params <- list(size = 0.3, alpha = 0.3, font.size = 3)
cao$sample.groups.palette <- c(Control = "#4575b4", AD = "#d73027")
```

The graph used is Seurat’s kNN graph on the Harmony embedding (`RNA_nn`,
symmetrized). The thin black lines on the maps below outline the regions
where the cells significant at an adjusted p of 0.05 are the majority
(`outline = TRUE`).

``` r
cowplot::plot_grid(cao$plotEmbedding(color.by = "cell.groups", title = "cell type"),
                   cao$plotEmbedding(color.by = "condition", title = "condition", mark.groups = FALSE, show.legend = TRUE), ncol = 2)
```

![](alz_cacoa_harmony_files/figure-gfm/embedding-1.png)<!-- -->

``` r
cao$setModel(~ condition, test = "condition")
```

    ## Model: ~condition      dispersion: ~condition
    ## Test: condition: AD vs Control  (reference 'Control': first level of the factor as set; to change: test = "condition: Control vs AD")
    ##        unadjusted; permutations: block (924 distinct); 12 of 12 samples can change their value; smallest p-value 0.0011
    ##        shift > 0: AD samples differ from Control samples in a common direction, beyond within-group variability
    ## Samples: 12 used.   Issues: none

## Cluster-free density shifts

The same map as in the main notebook: the default diffusion
neighborhoods (`scale = "auto"`, a scan over 1, 3 and 5 % of the cells),
each sample’s kernel density tested on the log scale with the model,
cells adjusted by the pooled permutation false-discovery rate.

``` r
cao$estimateCellDensity(method = "graph")
```

    ## kNN neighborhoods of at least 97 cells (8.1 per sample)

``` r
cao$estimateDiffCellDensity()
```

    ## diffusion neighborhoods at 1 %, 3 %, 5 % of the cells (auto: one scale is chosen)

    ## scale 0.01: heat kernel beta 18.6

    ## scale 0.03: heat kernel beta 25.2

    ## scale 0.05: heat kernel beta 29.4

    ## 7506 of 9758 cells at FDR 10 % (scale 5 % chosen)

``` r
cowplot::plot_grid(cao$plotEmbedding(color.by = "cell.groups", title = "cell type"),
                   cao$plotDiffCellDensity(what = "effect", title = "density change (log2 FC)", outline = TRUE),
                   cao$plotDiffCellDensity(title = "adjusted z", outline = TRUE), ncol = 3)
```

![](alz_cacoa_harmony_files/figure-gfm/density-1.png)<!-- -->

On the Harmony graph the scan again picks the 5 % scale (neighborhoods
of about 460 nuclei) and calls 7,506 of the 9,758 nuclei at FDR 10 % (77
%; 87 % on the Conos graph): every OPC (log2 fold change -1.2), 92 % of
the astrocytes (-1.1) and 79 % of the oligodendrocytes (+0.5), the three
types with the largest composition loadings, with per-type p-values
close to those on the Conos graph (0.008, adjusted 0.015; there
0.008-0.010, adjusted 0.019). The neurons, microglia and endothelial
cells are not significant in the per-type test on either graph. The
outlines enclose the three populations; on the oligodendrocytes they
leave out the part of the cluster where the change is smallest.

The kNN neighborhoods of the graph itself (`scale = NULL`: the stored
units of `estimateCellDensity()`, about 100 cells each) are the other
neighborhood definition available; they are one hop wide, so they cannot
leak the way a kernel can.

``` r
cao$estimateCellDensity(method = "graph", name = "cell.density.knn")
```

    ## kNN neighborhoods of at least 97 cells (8.1 per sample)

``` r
cao$estimateDiffCellDensity(name = "cell.density.knn", scale = NULL)
```

    ## 2958 of 9758 cells at FDR 10 %

``` r
cowplot::plot_grid(cao$plotDiffCellDensity(name = "cell.density.knn", what = "effect", title = "kNN units: density change (log2 FC)", outline = TRUE),
                   cao$plotDiffCellDensity(name = "cell.density.knn", title = "kNN units: adjusted z", outline = TRUE), ncol = 2)
```

![](alz_cacoa_harmony_files/figure-gfm/density-knn-1.png)<!-- -->

The kNN neighborhoods (97 nuclei) call fewer nuclei (2,958 at FDR 10 %):
73 % of the astrocytes and 76 % of the OPCs (log2 -1.8) but only 12 % of
the oligodendrocytes, whose moderate change (+0.5) needs the larger
neighborhoods to reach significance per nucleus; the per-type test calls
the astrocytes and OPCs (adjusted p 0.029), not the oligodendrocytes
(0.29).

## Cluster-free expression shifts

``` r
cao$estimateClusterFreeExpressionShifts(n.permutations = 499)
```

    ## 10373 genes (per neighborhood: those detected in >= 5 % of its cells)

    ## diffusion neighborhoods at 1 %, 3 %, 5 % of the cells (auto: one scale is chosen), at least 97 cells each (8.1 per sample)

    ## scale 0.01: heat kernel beta 18.6, 308 centers, neighborhoods of 102 cells (median, 98 smallest)

    ## scale 0.03: heat kernel beta 25.2, 125 centers, neighborhoods of 295 cells (median, 104 smallest)

    ## scale 0.05: heat kernel beta 29.4, 74 centers, neighborhoods of 491 cells (median, 152 smallest)

    ## 9709 of 9758 cells keep adjusted z >= 1.28 (scale 3 % chosen)

``` r
cao$plotClusterFreeExpressionShifts(adjusted = TRUE, font.size = c(2, 3), outline = TRUE)
```

![](alz_cacoa_harmony_files/figure-gfm/cluster-free-shifts-1.png)<!-- -->

The expression-shift map is as on the Conos graph: the scan picks 3 %,
97 % of the nuclei pass FDR 5 % here and 96 % there, and all six types
are significant here (five there: the endothelial cells, 72 nuclei, are
called at 0.017 instead of 0.062). The strongest shifts are in the OPCs
and astrocytes (median 1.6), then the endothelial cells (1.5), the
oligodendrocytes and microglia (1.4); the neurons have the smallest
(1.3) and only 47 % of their nuclei are called, where on the Conos graph
all of them are (shift 1.39).

``` r
cao$estimateClusterFreeSummary("cell.density")
cao$estimateClusterFreeSummary("cluster.free.expr.shifts")
cowplot::plot_grid(cao$plotClusterFreeSummary("cell.density"), cao$plotClusterFreeSummary("cluster.free.expr.shifts"), ncol = 2)
```

![](alz_cacoa_harmony_files/figure-gfm/cluster-free-summary-1.png)<!-- -->

Per type the two maps agree with the Conos notebook: three populations
change in density (OPCs, astrocytes, oligodendrocytes) and every type
changes in expression. On this small cohort of well-separated
populations the alignment makes little difference; the neurons are the
one population whose per-nucleus calls depend on the graph.

## Session info

``` r
sessionInfo()
```

    ## R version 4.2.2 Patched (2022-11-10 r83330)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.4    ggplot2_4.0.2  cacoa_0.5.0    Matrix_1.6-1.1
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] nlme_3.1-161           spatstat.sparse_3.1-0  matrixStats_1.5.0     
    ##   [4] RcppAnnoy_0.0.22       RColorBrewer_1.1-3     httr_1.4.7            
    ##   [7] sctransform_0.4.2      tools_4.2.2            R6_2.6.1              
    ##  [10] irlba_2.3.5.1          KernSmooth_2.23-20     colorspace_2.1-0      
    ##  [13] uwot_0.2.3             lazyeval_0.2.2         withr_3.0.2           
    ##  [16] sp_2.2-0               tidyselect_1.2.1       gridExtra_2.3         
    ##  [19] compiler_4.2.2         progressr_0.15.1       cli_3.6.5             
    ##  [22] spatstat.explore_3.5-2 plotly_4.11.0          labeling_0.4.3        
    ##  [25] Seurat_4.3.0.1         scales_1.4.0           spatstat.data_3.1-8   
    ##  [28] lmtest_0.9-40          S7_0.2.0               ggridges_0.5.7        
    ##  [31] pbapply_1.7-4          goftest_1.2-3          stringr_1.5.1         
    ##  [34] digest_0.6.37          spatstat.utils_3.1-5   rmarkdown_2.29        
    ##  [37] pkgconfig_2.0.3        htmltools_0.5.8.1      parallelly_1.45.1     
    ##  [40] fastmap_1.2.0          htmlwidgets_1.6.4      rlang_1.1.7           
    ##  [43] shiny_1.11.1           farver_2.1.2           generics_0.1.4        
    ##  [46] zoo_1.8-14             jsonlite_2.0.0         spatstat.random_3.4-1 
    ##  [49] ica_1.0-3              magrittr_2.0.3         patchwork_1.3.2       
    ##  [52] Rcpp_1.1.0             abind_1.4-8            ape_5.8-1             
    ##  [55] reticulate_1.43.0      lifecycle_1.0.4        stringi_1.8.7         
    ##  [58] yaml_2.3.10            MASS_7.3-58            Rtsne_0.17            
    ##  [61] plyr_1.8.9             grid_4.2.2             parallel_4.2.2        
    ##  [64] listenv_0.9.1          promises_1.3.3         ggrepel_0.9.6         
    ##  [67] deldir_2.0-4           miniUI_0.1.2           lattice_0.20-45       
    ##  [70] sccore_1.0.6           cowplot_1.2.0          splines_4.2.2         
    ##  [73] tensor_1.5.1           knitr_1.50             pillar_1.11.0         
    ##  [76] igraph_2.1.4           spatstat.geom_3.5-0    future.apply_1.20.0   
    ##  [79] reshape2_1.4.4         codetools_0.2-18       leiden_0.4.3          
    ##  [82] glue_1.8.0             evaluate_1.0.5         SeuratObject_4.1.3    
    ##  [85] spatstat.univar_3.1-4  data.table_1.17.8      png_0.1-8             
    ##  [88] vctrs_0.6.5            httpuv_1.6.16          polyclip_1.10-7       
    ##  [91] gtable_0.3.6           RANN_2.6.2             purrr_1.1.0           
    ##  [94] tidyr_1.3.1            scattermore_1.2        future_1.67.0         
    ##  [97] xfun_0.53              mime_0.13              xtable_1.8-4          
    ## [100] later_1.4.4            survival_3.4-0         viridisLite_0.4.2     
    ## [103] tibble_3.3.0           cluster_2.1.4          globals_0.18.0        
    ## [106] fitdistrplus_1.2-4     ROCR_1.0-11
