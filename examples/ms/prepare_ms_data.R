# Prepare the MS dataset (Schirmer et al. 2019) for ms_cacoa.Rmd from the objects of the Cacoa preprint
# (cacoaAnalysis/data/MS/con.rds: the Conos object with the preprint's filtering, 46,340 nuclei of 21 samples;
# cacoaAnalysis/data/MS/cao.rds: the annotation with the condition-specific EN-L2-3 and OL labels merged).
# Refreshes the R6 objects to the installed package versions and stores one list: con, cell.groups, sample.per.cell, embedding.
# usage: Rscript prepare_ms_data.R <out.rds> [<con.rds> <cao.rds>]
suppressMessages({library(conos); library(pagoda2); library(magrittr)})
args <- commandArgs(trailingOnly = TRUE); out.file <- args[1]
con.file <- if (length(args) > 1) args[2] else "/home/viktor_petukhov/projects/cacoaAnalysis/data/MS/con.rds"
cao.file <- if (length(args) > 2) args[3] else "/home/viktor_petukhov/projects/cacoaAnalysis/data/MS/cao.rds"
con <- readRDS(con.file)
for (n in names(con$samples)) con$samples[[n]] %<>% pagoda2::Pagoda2$new()
con <- Conos$new(con); con$graph <- igraph::upgrade_graph(con$graph)
ocao <- readRDS(cao.file); cell.groups <- droplevels(ocao$cell.groups); spc <- con$getDatasetPerCell()
cells <- igraph::V(con$graph)$name
cat(sprintf("%d samples; graph %d cells; annotation %d cells (%d levels), %d graph cells annotated; embedding %d\n", length(con$samples), length(cells), length(cell.groups), nlevels(cell.groups), sum(cells %in% names(cell.groups)), nrow(con$embedding)))
print(table(cell.groups)); print(table(spc[cells]))
saveRDS(list(con = con, cell.groups = cell.groups[intersect(names(cell.groups), cells)], sample.per.cell = droplevels(spc[cells]), embedding = con$embedding), out.file)
cat("wrote", out.file, "\n")
