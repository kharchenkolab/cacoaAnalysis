createPagoda <- function(cm, min.transcripts.per.cell=0, dub.scores=NULL, dub.threshold=0.25) {
  cm <- cm[-grep("^MT-", rownames(cm)), ]
  if (!is.null(dub.scores)) {
    cm <- cm[,dub.scores[colnames(cm)] < dub.threshold]
  }

  p2 <- pagoda2::Pagoda2$new(
    cm, n.cores=1, trim=10, log.scale=TRUE, min.cells.per.gene=0, verbose=FALSE,
    min.transcripts.per.cell=min.transcripts.per.cell
  )
  p2$adjustVariance(plot=FALSE, gam.k=10)
  p2$calculatePcaReduction(nPcs=100, n.odgenes=3000, maxit=1000)
  return(p2)
}

createConos <- function(p2s, sample.meta, cell.meta, n.cores=1, k=20, ...) {
  con <- conos::Conos$new(p2s, n.cores=n.cores)
  con$buildGraph(k=k, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN',
                 metric='angular', verbose=TRUE)
  con$embedGraph(method='UMAP', min.prob.lower=1e-4, verbose=FALSE, ...);

  con$misc$cell_metadata <- lapply(cell.meta, `[`, names(con$getDatasetPerCell()))
  con$misc$sample_metadata <- lapply(sample.meta, `[`, names(con$samples))
  return(con)
}

splitMatrixByFactor <- function(mat, fac) {
  tapply(1:ncol(mat), as.factor(fac[colnames(mat)]), function(ii) mat[,ii])
}

readOrCreate <- function(path, create.fun, force=FALSE) {
  if (!force && file.exists(path))
    return(readr::read_rds(path))

  res <- create.fun()
  readr::write_rds(res, path)
  return(res)
}
