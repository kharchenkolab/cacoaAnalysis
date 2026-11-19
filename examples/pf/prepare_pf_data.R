# Prepare the pulmonary fibrosis dataset (Habermann et al. 2020) for pf_cacoa.Rmd.
# Reads the Conos object of the Cacoa preprint (cacoaAnalysis/data/PF/con.rds: 22 pagoda2 samples, Control and IPF
# lungs only, one sample per lung, the paper's cell annotation in con$misc$cell_metadata$cellType and a UMAP embedding),
# refreshes the R6 objects to the installed pagoda2/conos versions, upgrades the stored igraph graph and stores the
# result as one Conos object.
# usage: Rscript prepare_pf_data.R <out.rds> [<cacoaAnalysis con.rds>]
suppressMessages({library(conos); library(pagoda2)})
args <- commandArgs(trailingOnly = TRUE)
out.file <- args[1]
in.file <- if (length(args) > 1) args[2] else "/home/viktor_petukhov/projects/cacoaAnalysis/data/PF/con.rds"

con <- readRDS(in.file)
for (n in names(con$samples)) con$samples[[n]] <- pagoda2::Pagoda2$new(con$samples[[n]])
con <- Conos$new(con)
con$graph <- igraph::upgrade_graph(con$graph)   # the stored graph predates the installed igraph
saveRDS(con, out.file)
