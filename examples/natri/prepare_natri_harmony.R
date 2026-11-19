# Harmony-aligned Seurat object of the Natri et al. 2024 pulmonary fibrosis lungs (GEO GSE227136) for natri_cacoa.Rmd: the cells, counts and
# annotation prepared in benchmarks/largedata (one library per donor, at most 2,000 cells per donor, genes expressed in >= 10 cells; 111 donors,
# 187,756 cells) on a kNN graph of a Harmony embedding (see ../prepare_harmony_common.R). The sample metadata is saved next to it.
# usage: Rscript prepare_natri_harmony.R <inputs_natri.rds> <sample_meta_natri.rds> <out.rds> [n.cores]
suppressPackageStartupMessages(library(Matrix))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 3) as.integer(args[4]) else 8
x <- readRDS(args[1]); meta <- readRDS(args[2]); stopifnot(all(levels(droplevels(x$sample.per.cell)) %in% rownames(meta)))
saveRDS(meta, sub("\\.rds$", "_sample_meta.rds", args[3]))
harmonySeurat(x$counts, x$sample.per.cell, x$cell.groups, args[3], n.cores = n.cores)
