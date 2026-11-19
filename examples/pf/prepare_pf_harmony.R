# Prepare a Harmony-aligned Seurat object of the pulmonary fibrosis dataset for pf_cacoa_harmony.Rmd.
# Takes the same cells, counts and annotation as the Conos object of pf_cacoa.Rmd (prepare_pf_data.R), and replaces the
# Conos alignment graph by a kNN graph on a Harmony embedding: log-normalization, 2,000 variable genes, 50 PCs, Harmony by
# sample, 30 nearest neighbors on 30 Harmony dimensions (Seurat's RNA_nn / RNA_snn graphs) and a UMAP on the same dimensions.
# usage: Rscript prepare_pf_harmony.R <pf_con_refreshed.rds> <out.rds> [n.cores]
suppressMessages({library(conos); library(Seurat); library(magrittr)})
args <- commandArgs(trailingOnly = TRUE)
in.file <- args[1]; out.file <- args[2]; n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
set.seed(1)
con <- readRDS(in.file)
cm <- con$getJointCountMatrix(raw = TRUE)                       # cells x genes, the raw counts of all 22 samples
spc <- con$getDatasetPerCell(); spc <- setNames(as.character(spc[rownames(cm)]), rownames(cm))
cmeta <- con$misc$cell_metadata; grp <- setNames(as.character(cmeta$cellType), cmeta$cell)[rownames(cm)]
rm(con); invisible(gc())
seu <- CreateSeuratObject(counts = Matrix::t(cm), min.cells = 0, min.features = 0)
seu$sample <- spc[colnames(seu)]; seu$celltype <- grp[colnames(seu)]
seu <- NormalizeData(seu, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:50, reduction.save = "harmony", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, k.param = 30, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
saveRDS(seu, out.file)
message(sprintf("saved %s: %d cells, %d genes", out.file, ncol(seu), nrow(seu)))
