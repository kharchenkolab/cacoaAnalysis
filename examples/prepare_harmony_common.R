# Shared by prepare_<dataset>_harmony.R: a Harmony-aligned Seurat object from a cells x genes raw count matrix, the sample of
# each cell and the cell annotation. Log-normalization, 2,000 variable genes, 50 PCs, Harmony by sample, 30 nearest neighbors
# on 30 Harmony dimensions (Seurat's RNA_nn / RNA_snn graphs) and a UMAP on the same dimensions.
harmonySeurat <- function(cm, sample.per.cell, cell.groups, out.file, n.cores = 8, seed = 1) {
  suppressMessages({library(Seurat); library(magrittr)})
  set.seed(seed)
  spc <- as.character(sample.per.cell[rownames(cm)]); grp <- as.character(cell.groups[rownames(cm)])
  ok <- !is.na(spc) & !is.na(grp)                                   # cells of the count matrices that are not in the analysis (no sample or annotation)
  if (any(!ok)) message(sprintf("dropping %d of %d cells without a sample or an annotation", sum(!ok), length(ok)))
  cm <- cm[ok, , drop = FALSE]; spc <- spc[ok]; grp <- grp[ok]; cells <- rownames(cm)
  message(sprintf("counts: %d cells x %d genes", nrow(cm), ncol(cm)))
  seu <- CreateSeuratObject(counts = Matrix::t(cm), min.cells = 0, min.features = 0)
  seu$sample <- spc; seu$celltype <- grp
  seu <- NormalizeData(seu, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
  message("PCA done")
  seu <- harmony::RunHarmony(seu, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:50, reduction.save = "harmony", verbose = FALSE)
  message("Harmony done")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, k.param = 30, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
  saveRDS(seu, out.file)
  message(sprintf("saved %s: %d cells, %d genes, %d samples", out.file, ncol(seu), nrow(seu), length(unique(seu$sample))))
  invisible(seu)
}
