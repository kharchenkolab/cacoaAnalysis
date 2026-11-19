# Harmony-aligned Seurat object of the MS dataset for ms_cacoa_harmony.Rmd: the cells, counts and annotation of the Conos
# object of ms_cacoa.Rmd (ms_refreshed.rds) on a kNN graph of a Harmony embedding (see ../prepare_harmony_common.R).
# usage: Rscript prepare_ms_harmony.R <ms_refreshed.rds> <out.rds> [n.cores]
suppressMessages(library(conos)); source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
ms <- readRDS(args[1]); con <- ms$con
cms <- lapply(con$samples, function(s) s$misc$rawCounts)      # cells x genes per sample; the gene sets differ between samples
genes <- Reduce(intersect, lapply(cms, colnames))
cm <- do.call(rbind, lapply(cms, function(m) m[, genes, drop = FALSE]))
spc <- con$getDatasetPerCell()
harmonySeurat(cm, spc, ms$cell.groups, args[2], n.cores = n.cores)
