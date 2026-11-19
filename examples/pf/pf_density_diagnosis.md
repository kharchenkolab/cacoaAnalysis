# PF: why the cluster-free density map calls no cell while the composition changes clearly

Question (2026-10-29): the density map of the PF notebook reports "0 of 89,326 cells at FDR 10 % (scale 1 % chosen)" although the
composition analysis calls the basal cells (0.2 % of the cells in the controls, 3.4 % in IPF), `SCGB3A2+` (0.5 to 5.6 %),
`KRT5-/KRT17+` (0 to 0.3 %) and AT1 cells (1.3 to 0.1 %) at adjusted p 0.02-0.03. Scripts and logs: scratchpad `pfdens/`
(`check1-3.R`, `modprior.R`, `diffscale.R`, `kernel.R`, `d0max.R`, `knn_vs_diffusion_{ms,scc}.R`); stored result `ex_pf/results_v27.rds`.

## What is not wrong

- The samples are aligned: the columns of the density matrix are the rows of the model's design, in the same order (10 Control, 12 IPF).
- The map points the right way and ranks the right types. Per cell type, the group-level permutation test of the stored map gives
  basal 0.006, mesothelial 0.006, `MUC5AC+ High` 0.008, `MUC5B+` 0.008, `KRT5-/KRT17+` 0.024 (31 types, BH 0.06-0.11); the
  composition's raw p-values for its calls are in the same range (0.002-0.004 before the adjustment over 31 types). The loss is at
  the cell level, where the pooled FDR runs over 89,326 cells.

## Where the per-cell signal goes

1. **The heat kernel leaks out of the small cell types on this Conos graph** (89,326 cells, 3.0 M edges, median degree 39). Kernel
   columns of 25 probe cells per type (`kernel.R`): at the 1 % scale (nominal 893 cells, effective size 854) a basal probe keeps 65 %
   of its kernel mass on basal cells, an AT1 probe 21 %, a `KRT5-/KRT17+` probe 16 % (macrophages 95 %); 90 % of the mass is spread
   over 5,700-8,000 cells, 6-9 times the nominal size. At 0.3 % (268 cells) AT1 keeps 27 %. On the MS graph the same probes keep
   80-100 % of their mass on their own type at 1 %. The per-cell effects are therefore regional averages: AT1 cells read +0.2 log2
   instead of -0.9, every fibroblast subtype reads about +0.8 like the expanding `PLIN2+` fibroblasts next to them, and the basal
   fold change shrinks from 4.2 log2 (composition) to 2.3 (the 220-cell kNN neighborhood: 14 % of its cells are control cells,
   matched in by the alignment) to 1.2 (diffusion at 1 %).
2. **The trend moderation degenerated.** limma's `fitFDist` on the label-free variances returned an infinite prior df (capped at
   1e6), so every cell's own variance was replaced by the trend value, which is larger than the within-arm variance of exactly the
   cells with a difference (basal 2.7 vs 4.4; t shrinks by 0.78). On the kNN units in a direct call: 1,027 cells at FDR 10 % with
   the moderation, 3,510 without; the pooled null without it matches the parametric t (pooled p 0.0024 vs 0.0028 for basal cells).

## Results of the alternatives (199 permutations, `cao$estimateDiffCellDensity()` through the class)

| neighborhoods | cells at FDR 10 % | basal | AT1 | `KRT5-/KRT17+` | `SCGB3A2+` |
|---|---|---|---|---|---|
| diffusion 1 % (default choice) | 0 (min q 0.28) | 0 | 0 | 0 | 0 |
| diffusion 0.3 % | 0 (min q 0.23) | 0 | 0 | 0 | 0 |
| diffusion 1 %, prior df capped at 4 | 0 (min q 0.18) | 0 | 0 | 0 | 0 |
| kNN units (`scale = NULL`, 220 cells) | 4,599 | 82 % | 20 % | 27 % | 26 % |

Direct calls on the kNN units with the own prior fitter and a cap on the prior df: 964 cells uncapped, 2,160 (cap 40), 2,522 (20),
3,007 (10), 3,301 (4), 3,361 (2), 3,510 unmoderated. The class's kNN path adds the neighbor averaging of the statistic (4,599).

On the other cohorts the kernel is the better unit: MS 19,467 cells at FDR 10 % on the diffusion scan (3 % chosen) against 9,448 on
the kNN units; SCC 17,735 (5 %) against 11,052.

## Why the kernel leaks on PF (graph structure, `graphstruct.R`)

Both Conos graphs consist of cross-sample edges (mutual nearest neighbors between samples; within-sample edges carry 1-2 % of
the weight; edge weights 0.95 median on both). A cell of a population that a sample lacks is matched to that sample's closest
other population, and those edges weigh as much as true matches. Per cell type, the fraction of edge weight leaving the type, and
the share of those cross-type edges whose other end lies in a sample with fewer than 10 cells of the type:

| PF type | cells | samples with it (of 22) | weight leaving the type | of which to samples lacking it |
|---|---|---|---|---|
| `KRT5-/KRT17+` | 387 | 7 | 0.62 | 0.65 |
| `PLIN2+` fibroblasts | 1,173 | 4 | 0.58 | 0.90 |
| `SCGB3A2+` | 2,921 | 14 | 0.49 | 0.41 |
| AT1 | 568 | 10 | 0.40 | 0.52 |
| mesothelial | 172 | 5 | 0.36 | 0.88 |
| basal | 1,851 | 14 | 0.26 | 0.26 |
| macrophages | 28,511 | 22 | 0.10 | 0.00 |
| ciliated | 13,976 | 21 | 0.02 | 0.01 |

Eleven of the 31 PF types lose more than half their edge weight to other types (weighted mean 0.18). On MS the types that carry
the result are present in 20-21 of 21 samples and lose 1-27 % of their weight, almost none of it to samples lacking them
(EN-L2-3 0.24 / 0.007, microglia 0.06 / 0.53 of a small amount, OL 0.01, astrocytes 0.01); only the tiny T, B, phagocyte and
stromal populations (43-132 cells, 1-4 samples) leak like the PF populations. The heat kernel compounds the leak: once the heat
has crossed into another type it spreads through that type's dense within-type edges, so a kernel of 268 cells keeps 27 % of its
mass on AT1 cells while the one-hop kNN unit of 220 cells stays mostly within the type. It is not the connectivity as such
(median degree 39 on PF, 77 on MS) but the alignment of populations that exist in only part of the samples: the IPF-emergent
epithelium and fibroblasts, exactly the populations the composition calls.

## The same cells on a Harmony-aligned graph (`prepare_pf_harmony.R`, `pf_cacoa_harmony.Rmd`)

A kNN graph (k = 30) on a Harmony embedding of the same counts and annotation (2,000 variable genes, 50 PCs, Harmony by
sample, 30 dimensions): 89,326 cells, 2.1 M edges, median degree 42, one component, and 31 % of the edge weight within
samples (2 % on the Conos graph). The kernel stays home: at 1 % a basal probe keeps 92 % of its mass on basal cells (Conos
65 %), AT1 72 % (21 %), `KRT5-/KRT17+` 51 % (16 %), and 90 % of the mass lies within 520-950 cells (Conos 5,700-8,000). The
edge weight leaving the type drops to 0.29 for `KRT5-/KRT17+` (from 0.62), 0.27 for the `PLIN2+` fibroblasts (0.58), 0.29
for AT1 (0.40) and 0.09 for basal (0.26). The density map (199 permutations, through the class) calls 11,700 cells at FDR
10 % on the diffusion scan (1 % chosen; basal 96 %, `KRT5-/KRT17+` 80 %, `SCGB3A2+` 60 %, AT1 8 %) and 26,836 on the kNN
units (basal 96 %, `KRT5-/KRT17+` 89 %, AT1 53 %, AT2 26 %); the effects are less compressed (basal log2 1.85 on the kernel,
2.79 on the kNN units). So the null map is a property of the Conos alignment of this cohort, not of the statistic.

## What was changed in the package (rlib_ss36/37, `cacoa/NEWS.md` "Statistics")

- `fitFPrior()` fits the scaled-F prior itself (moments of the log variances, lowess trend on the mean; agrees with limma's
  `fitFDist` on simulated variances) and caps the prior df; `estimateDiffCellDensity(moderate.d0.max = 4)` is the default, so a
  cell's own variance keeps at least 83 % of the weight at 20 residual df. limma is no longer needed for the moderation (it stays a
  suggested package for limma-voom DE).

## Open (decision)

- The default scale scan over 1/3/5 % cannot rescue PF: the kernel does not localize on this graph at any scale. Including the
  kNN units as one more candidate of the auto scan, chosen by the existing rule (the unit calling the most cells at FDR 10 %),
  would pick the kNN units on PF (4,599 vs 0) and leave MS and SCC on the kernel. A kernel-purity diagnostic (mass of a probe's
  kernel on its own annotation) would flag PF without the annotation entering the test.
- Until then the PF notebook's density section describes a power problem; with `scale = NULL` it would show the basal, `SCGB3A2+`,
  `KRT5-/KRT17+` and AT1 changes.
