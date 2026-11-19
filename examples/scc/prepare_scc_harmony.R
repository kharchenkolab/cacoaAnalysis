# Harmony-aligned Seurat object of the SCC dataset for scc_cacoa_harmony.Rmd: the cells, counts and annotation of the Conos
# object of scc_cacoa.Rmd (scc_refreshed.rds) on a kNN graph of a Harmony embedding (see ../prepare_harmony_common.R).
# usage: Rscript prepare_scc_harmony.R <scc_refreshed.rds> <out.rds> [n.cores]
suppressMessages(library(conos)); source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
scc <- readRDS(args[1]); con <- scc$con
cm <- con$getJointCountMatrix(raw = TRUE); spc <- scc$sample.per.cell
harmonySeurat(cm, spc, scc$cell.groups, args[2], n.cores = n.cores)
