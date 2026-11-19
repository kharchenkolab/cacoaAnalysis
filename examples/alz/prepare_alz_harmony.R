# Harmony-aligned Seurat object of the ALZ dataset for alz_cacoa_harmony.Rmd: the cells, counts and annotation of the Cacoa
# object of the previous package version used by alz_cacoa.Rmd, on a kNN graph of a Harmony embedding (see ../prepare_harmony_common.R).
# usage: Rscript prepare_alz_harmony.R <cao.rds> <out.rds> [n.cores]
suppressMessages({library(conos); library(pagoda2)}); source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
ocao <- readRDS(args[1])
for (n in names(ocao$data.object$samples)) ocao$data.object$samples[[n]] <- pagoda2::Pagoda2$new(ocao$data.object$samples[[n]])
cms <- lapply(ocao$data.object$samples, function(s) s$misc$rawCounts)      # cells x genes per sample; the gene sets differ between samples
genes <- Reduce(intersect, lapply(cms, colnames))
cm <- do.call(rbind, lapply(cms, function(m) m[, genes, drop = FALSE]))
harmonySeurat(cm, ocao$sample.per.cell, ocao$cell.groups, args[2], n.cores = n.cores)
