# Harmony-aligned Seurat object of a subset of OneK1K (Yazar et al. 2022, CELLxGENE; 981 healthy donors aged 19-97) for onek1k_cacoa.Rmd: 200
# donors drawn evenly over age (ten age bins of equal donor count) and sex (seed 1), at most 800 cells per donor (uniformly within each donor),
# genes expressed in >= 10 cells, the authors' cell types with the merging of benchmarks/largedata (rare types folded into their neighbours,
# doublets and erythrocytes removed), on a kNN graph of a Harmony embedding (see ../prepare_harmony_common.R). The donors' metadata
# (pool, sex, age) is saved next to it.
# usage: Rscript prepare_onek1k_harmony.R <counts_onek1k.h5> <meta_onek1k.rds> <selected_cells.tsv> <out.rds> [n.cores] [donors = 200] [cells per donor = 800]
suppressPackageStartupMessages({library(rhdf5); library(Matrix)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "prepare_harmony_common.R"))
args <- commandArgs(trailingOnly = TRUE); n.cores <- if (length(args) > 4) as.integer(args[5]) else 8
n.donors <- if (length(args) > 5) as.integer(args[6]) else 200L; cap <- if (length(args) > 6) as.integer(args[7]) else 800L
mt <- readRDS(args[2]); meta <- mt$sample.meta; sel <- read.delim(args[3]); stopifnot(all(c("cell", "donor") %in% names(sel)))
# donors: ten age bins of equal donor count, then an equal share of each sex within a bin where possible
set.seed(1); meta$age.bin <- cut(rank(meta$age, ties.method = "first"), 10, labels = FALSE)
per.bin <- n.donors / 10
pick <- unlist(lapply(split(rownames(meta), meta$age.bin), function(d) { s <- split(d, meta[d, "sex"]); k <- floor(per.bin / 2)
  a <- unlist(lapply(s, function(v) sample(v, min(length(v), k)))); rest <- setdiff(d, a); c(a, sample(rest, per.bin - length(a))) }))
meta <- meta[pick, setdiff(names(meta), "age.bin")]; message(sprintf("%d donors; age %d-%d; %s", nrow(meta), min(meta$age), max(meta$age), paste(names(table(meta$sex)), table(meta$sex), collapse = " / ")))
h <- args[1]; shape <- as.integer(h5read(h, "shape")); cells <- as.character(h5read(h, "cells")); sym <- as.character(h5read(h, "gene_symbol")); gid <- as.character(h5read(h, "gene_id"))
spc.all <- setNames(as.character(sel$donor), sel$cell)
take <- which(cells %in% names(spc.all) & spc.all[cells] %in% rownames(meta)); take <- take[!is.na(take)]
by.donor <- split(take, spc.all[cells[take]]); take <- sort(unlist(lapply(by.donor, function(i) if (length(i) > cap) sample(i, cap) else i), use.names = FALSE))
indptr <- as.numeric(h5read(h, "indptr")); nn <- diff(indptr)[take]
idx <- unlist(lapply(take, function(r) if (indptr[r + 1] > indptr[r]) seq(indptr[r] + 1, indptr[r + 1]) else NULL))   # the selected cells' entries (rows of a cells x genes CSR matrix)
ind <- h5read(h, "indices")[idx] + 1L; val <- as.numeric(h5read(h, "data")[idx]); invisible(gc())
cm <- sparseMatrix(i = rep.int(seq_along(take), nn), j = ind, x = val, dims = c(length(take), shape[2]))
sym[duplicated(sym) | sym == ""] <- gid[duplicated(sym) | sym == ""]; dimnames(cm) <- list(cells[take], sym)
ct <- mt$cell.type[rownames(cm)]; merge <- c(ASDC = "cDC2", cDC1 = "cDC2", Doublet = NA, Eryth = NA, `CD8 Proliferating` = "CD8 TEM", `CD4 Proliferating` = "CD4 TEM", `NK Proliferating` = "NK", ILC = "NK_CD56bright", dnT = "gdT")
ct <- ifelse(ct %in% names(merge), merge[ct], ct); names(ct) <- rownames(cm)
cm <- cm[!is.na(ct), , drop = FALSE]; cm <- cm[, Matrix::colSums(cm > 0) >= 10, drop = FALSE]
message(sprintf("%d cells of %d donors, %d genes", nrow(cm), length(unique(spc.all[rownames(cm)])), ncol(cm)))
saveRDS(meta, sub("\\.rds$", "_sample_meta.rds", args[4]))
harmonySeurat(cm, spc.all[rownames(cm)], ct, args[4], n.cores = n.cores)
