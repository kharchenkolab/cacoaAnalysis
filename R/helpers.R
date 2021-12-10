#' @import magrittr
#' @import ggplot2

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
  con$embedGraph(method='UMAP', min.prob.lower=1e-5, verbose=FALSE, ...);

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

initializeCacoa <- function(dataset.name, ...){
  dataset.name %<>% toupper()
  con <- dataorganizer::DataPath(dataset.name, "con.rds") %>%
    readr::read_rds() %>% conos::Conos$new()
  ref.level <- 'Control'

  if (dataset.name == 'AZ') {
    target.level <- 'Alzheimer'
    sample.groups <- con$misc$sample_metadata$sample %>%
      {setNames(ifelse(startsWith(., 'A'), target.level, ref.level), .)}
  } else if (dataset.name == 'ASD') {
    target.level <- 'ASD'
    sample.groups <- con$misc$sample_metadata$diagnosis
  } else if (dataset.name == 'EP') {
    target.level <- 'Epilepsy'
    sample.groups <- con$misc$sample_metadata$Alias %>%
      {setNames(ifelse(startsWith(., 'E'), target.level, ref.level), .)}
  } else if (dataset.name == 'MS') {
    target.level <- 'MS'
    sample.groups <- con$misc$sample_metadata$diagnosis
  } else if (dataset.name == 'PF') {
    target.level <- 'IPF'
    sample.groups <- con$misc$sample_metadata$Diagnosis
  } else if (dataset.name == 'SCC') {
    target.level <- 'Tumor'
    sample.groups <- names(con$samples) %>%
      {setNames(ifelse(grepl('Normal', .), ref.level, target.level), .)}
  }

  cao <- cacoa::Cacoa$new(
    con, sample.groups=sample.groups, cell.groups=con$misc$cell_metadata$cellType,
    target.level=target.level, ref.level=ref.level, ...
  )

  bs <- if (length(levels(cao$cell.groups)) > 20) 10 else 16
  cao$plot.theme <- theme_bw(base_size = bs) +
    theme(plot.title=element_text(hjust = 0.5),
          legend.background=element_rect(fill=alpha("white", 0.2)),
          legend.text=element_text(size=12),
          legend.margin=margin(6, 6, 4, 1, 'pt'),
          plot.margin=margin())

  set.seed(239)
  cao$cell.groups.palette <- levels(cao$cell.groups) %>%
    {setNames(sample(brewerPalette("Paired")(length(.))), .)}
  cao$sample.groups.palette <- c("#d73027", "#4575b4") %>%
    setNames(c(cao$target.level, cao$ref.level))

  return(cao)
}

runCacoaAnalyses <- function(cao, cluster.free.de=FALSE, cluster.based.de=FALSE, cell.density=c('graph', 'kde'),
                             verbose=2, n.cores.ontology=1) {
  verb1 <- (verbose > 0)
  verb2 <- (verbose > 1)
  cao$estimateExpressionShiftMagnitudes(n.permutations=2500, verbose=verb2)
  cao$estimateCellLoadings()

  for (met in c('kde', 'graph')) {
    if (!(met %in% cell.density))
      next

    cn <- paste0('cell.density.', met)
    cao$estimateCellDensity(method=met, estimate.variation=FALSE, verbose=verb2, name=cn, beta=10)
    cao$estimateDiffCellDensity(type='wilcox', adjust.pvalues=TRUE, verbose=verb2, n.permutations=500, name=cn)
    cao$estimateDiffCellDensity(type='subtract', adjust.pvalues=FALSE, verbose=verb2, name=cn)
  }

  cao$estimateClusterFreeExpressionShifts(n.top.genes=3000, gene.selection="expression", verbose=verb1)

  if (cluster.free.de) {
    cao$estimateClusterFreeDE(n.top.genes=1500, min.expr.frac=0.01, adjust.pvalues=TRUE, smooth=TRUE, verbose=verb1)
    cao$smoothClusterFreeZScores(progress.chunks=10, z.adj=TRUE, verbose=verb1)
  }

  if (cluster.based.de) {
    cao$estimateDEPerCellType(independent.filtering=TRUE, test='DESeq2.Wald', verbose=verb1)
    cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=verb1, n.cores=n.cores.ontology)
  }

  return(cao)
}

figurePath <- function(...) dataorganizer::OutputPath("figures", ...)
