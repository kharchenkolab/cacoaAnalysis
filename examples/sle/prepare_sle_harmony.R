# Harmony-aligned Seurat object of the Perez et al. 2022 lupus PBMC data (CELLxGENE, 261 donors: 162 SLE, 99 healthy) for sle_cacoa.Rmd: the cells,
# counts and annotation prepared in benchmarks/largedata (one library group per donor, at most 2,000 cells per donor), subsampled to at most
# 600 cells per donor (uniformly within each donor, seed 1) to keep the notebook tractable, on a kNN graph of a Harmony embedding
# (see ../prepare_harmony_common.R). The sample metadata is saved next to it.
# usage: Rscript prepare_sle_harmony.R <inputs_sle.rds> <out.rds> [n.cores] [cells per donor = 600]
suppressPackageStartupMessages(library(Matrix))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 2) as.integer(args[3]) else 8; cap <- if (length(args) > 3) as.integer(args[4]) else 600L
x <- readRDS(args[1]); spc <- x$sample.per.cell; set.seed(1)
keep <- unlist(lapply(split(names(spc), spc), function(cl) if (length(cl) > cap) sample(cl, cap) else cl), use.names = FALSE)
message(sprintf("%d of %d cells kept (at most %d per donor)", length(keep), length(spc), cap))
cm <- x$counts[keep, , drop = FALSE]; cm <- cm[, Matrix::colSums(cm > 0) >= 10, drop = FALSE]
saveRDS(x$sample.meta, sub("\\.rds$", "_sample_meta.rds", args[2]))
harmonySeurat(cm, spc, x$cell.groups, args[2], n.cores = n.cores)
