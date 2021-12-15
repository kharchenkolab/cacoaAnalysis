#' @import parallel magrittr dplyr Matrix
NULL

## Prepare data

prepareSCEForSimulaton <- function(cm, sample.per.cell, ref.group, cell.groups, sample.groups,
                                   subgroups=unique(cell.groups)) {
  meta.df <- cell.groups %>% .[. %in% subgroups] %>%
    data.frame(cluster_id=as.factor(.)) %>%
    cbind(sample_id=sample.per.cell[rownames(.)]) %>%
    cbind(group_id=as.factor(sample.groups[as.character(.$sample_id)]))

  sce <- list(counts=cm[,rownames(meta.df)]) %>%
    SingleCellExperiment::SingleCellExperiment(colData=meta.df)

  prep <- muscat::prepSim(sce, group_keep=ref.group)
  offset <- prep@colData$offset
  prep@colData <- sce@colData[rownames(prep@colData),]
  prep@colData$offset <- offset

  return(list(raw=sce, prep=prep))
}

cacoaFromSim <- function(sim, n.cores, ref.level="A", target.level=setdiff(unique(sim$sample.groups), ref.level)) {
  sim %$% cacoa::Cacoa$new(cm, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell,
                           ref.level=ref.level, target.level=target.level, n.cores=n.cores, embedding=NULL)
}

## Generate simulations

convertSCEToSim <- function(sce, params, cluster.id=paste(params, collapse="_"), store.sce=FALSE) {
  cm <- sce@assays@data$counts %>% as("dgCMatrix")
  colnames(cm) %<>% paste0("_", cluster.id)
  meta <- sce@colData %>% as.data.frame() %>% set_rownames(paste0(rownames(.), "_", cluster.id))
  meta$sample_id %<>% droplevels()

  sample.groups <- meta %$% split(group_id, sample_id) %>% sapply(unique)
  cell.groups <- rep(paste0(cluster.id), nrow(meta)) %>%
    setNames(rownames(meta))
  sample.per.cell <- setNames(meta$sample_id, rownames(meta))

  params$cluster.id <- cluster.id
  params <- params[!sapply(params, is.null)]
  res <- list(
    cm=cm, sample.groups=sample.groups, cell.groups=cell.groups,
    sample.per.cell=sample.per.cell, params=as_tibble(params),
    gene.info=sce@metadata$gene_info
  )

  if (store.sce) {
    res$sce <- sce
  }

  return(res)
}

subsampleSCE <- function(sce, sample.frac=1.0, cell.frac=1.0, sample.condition=NULL) {
  params <- list(sample.frac=sample.frac, cell.frac=cell.frac) %>% lapply(signif, 3)
  params$sample.condition <- sample.condition

  if (cell.frac < 1.0) {
    subs.ids <- split(1:ncol(sce), sce$sample_id) %>%
      lapply(function(ids) sample(ids, cell.frac * length(ids))) %>%
      do.call(c, .)
    sce <- sce[, subs.ids]
  }

  if (sample.frac < 1) {
    sample.per.group <- sce@colData %>% as.data.frame() %$%
      split(group_id, sample_id) %>% sapply(unique) %>% {split(names(.), .)}
    if (is.null(sample.condition)) {
      sample.per.group %<>% lapply(function(ids) sample(ids, round(sample.frac * length(ids))))
    } else {
      sample.per.group[[sample.condition]] %<>% sample(round(sample.frac * length(.)))
    }

    sce <- sce[, sce$sample_id %in% do.call(c, sample.per.group)]
  }

  return(convertSCEToSim(sce, params))
}

getDEFracVector <- function(de.fraction, de.type){
  de.types <- c('ep', 'de', 'dp', 'dm', 'db')
  de.vec <- rep(0, 6)
  de.vec[1] <- 1 - de.fraction
  de.vec[which(de.types == de.type) + 1] <- de.fraction
  return(de.vec)
}

generateSimData <- function(sce, ns=length(unique(sce$sample_id)), nc=300, de.frac=0.05, de.type="de", lfc=1, rep=1,
                            cluster.id=NULL, nc.fixed=FALSE, ns.condition=NULL, store.sce=FALSE, suffix=NULL, ...) {
  params <- list(ns=ns, nc=nc, de.frac=de.frac, de.type=de.type, lfc=lfc, rep=rep, suffix=suffix)

  if (is.null(cluster.id)) cluster.id <- paste(params, collapse="_")

  if (!is.null(ns.condition)) {
    ns.cond <- ns
    ns <- length(unique(sce$sample_id))
  }

  if (!nc.fixed) nc <- nc * ns * 2

  de.vec <- getDEFracVector(de.frac, de.type=de.type)
  sim <- muscat::simData(sce, nk=1, ns=ns, nc=nc, p_dd=de.vec, lfc=lfc, ...)

  if (!is.null(ns.condition)) {
    samp.ids <- sim$sample_id[sim$group_id == ns.condition] %>% unique() %>% sample(ns.cond) %>%
      c(unique(sim$sample_id[sim$group_id != ns.condition]))
    sim <- sim[, sim$sample_id %in% samp.ids]
  }

  return(convertSCEToSim(sim, params, cluster.id=cluster.id, store.sce=store.sce))
}

joinSims <- function(sims) {
  sims <- unname(sims)
  params <- lapply(sims, `[[`, "params") %>% do.call(rbind, .)
  gene.info <- lapply(sims, function(sm) {
    if ("cluster.id" %in% colnames(sm$gene.info)) sm$gene.info else cbind(sm$gene.info, cluster.id=sm$params$cluster.id)
  }) %>% do.call(rbind, .)
  sces <- NULL
  if (!is.null(sims[[1]]$sce)) {
    sces <- lapply(sims, `[[`, "sce")
    sims[[1]]$sce <- NULL
  }

  res <- names(sims[[1]]) %>% setdiff(c("params", "gene.info")) %>%
    setNames(., .) %>% lapply(function(n) do.call(c, lapply(sims, `[[`, n)))
  res$cm %<>% sccore::mergeCountMatrices()
  res$sample.groups %<>% split(names(.)) %>% sapply(unique)
  res$params <- params
  res$gene.info <- gene.info
  res$sce <- sces

  return(res)
}

generateSims <- function(sce, lfc=1, n.samples=length(unique(sce$sample_id)), n.cells=300, de.frac=0.05, de.type="de", n.cores=1, n.repeats=1, verbose=TRUE, ...) {
  nc2 <- max(n.cores %/% length(lfc), 1)
  nc3 <- max(nc2 %/% length(n.samples), 1)
  nc4 <- max(nc3 %/% length(n.cells), 1)
  nc5 <- max(nc4 %/% length(de.frac), 1)
  nc6 <- max(nc5 %/% length(de.type), 1)

  if (verbose) {
    cat(n.cores, nc2, nc3, nc4, nc5, nc6)
  }

  res <- mclapply(lfc, function(lf) {
    mclapply(n.samples, function(ns) {
      mclapply(n.cells, function(nc) {
        mclapply(de.frac, function(def) {
          mclapply(de.type, function(det) {
            mclapply(1:n.repeats, function(rep) {
              generateSimData(sce, ns=ns, nc=nc, lfc=lf, de.frac=def, de.type=det, rep=rep, ...)
            }, mc.preschedule=TRUE, mc.cores=nc6, mc.allow.recursive=TRUE)
          }, mc.preschedule=TRUE, mc.cores=nc5, mc.allow.recursive=TRUE) %>% do.call(c, .)
        }, mc.preschedule=TRUE, mc.cores=nc4, mc.allow.recursive=TRUE) %>% do.call(c, .)
      }, mc.preschedule=TRUE, mc.cores=nc3, mc.allow.recursive=TRUE) %>% do.call(c, .)
    }, mc.preschedule=TRUE, mc.cores=nc2, mc.allow.recursive=TRUE) %>% do.call(c, .)
  }, mc.preschedule=TRUE, mc.cores=n.cores, mc.allow.recursive=TRUE) %>% do.call(c, .)
  if (length(res) == 1) return(res[[1]])
  return(joinSims(res))
}

## Plotting

prepareExpressionDistDf <- function(dist.per.type, params=NULL, clust.order=NULL) {
  if (is.null(clust.order) && !is.null(params)) {
    clust.order <- arrange(params, lfc, de.type, de.frac, ns, nc)$cluster.id
  }

  res <- names(dist.per.type) %>%
    lapply(function(n) tibble(Type=n, value=dist.per.type[[n]])) %>%
    do.call(rbind, .) %>% na.omit() %>% mutate(Type=as.factor(Type))

  if (!is.null(clust.order)) res %<>% mutate(Type=droplevels(factor(Type, levels=clust.order)))
  return(res)
}

## Cacoa metrics

getExprShift <- function(cao.obj, dist='cor', n.subsamples=50, params=NULL, ...) {
  res <- cao.obj$estimateExpressionShiftMagnitudes(dist=dist, n.subsamples=n.subsamples, ...)
  dist.per.type <- res$dist.df %$% split(value, Type)
  if (is.null(params)) return(dist.per.type)
  return(prepareExpressionDistDf(dist.per.type, params=params))
}

getCommonExprShift <- function(cao.obj, params=NULL, ...){
  res <- cao.obj$estimateCommonExpressionShiftMagnitudes(...)
  cell.types <- names(res[[1]])
  dist.per.type <- lapply(setNames(cell.types, cell.types), function(n) {
    lapply(res, `[[`, n) %>% do.call(rbind, .) %>% colMeans()
  })
  if (is.null(params)) return(dist.per.type)
  return(prepareExpressionDistDf(dist.per.type, params=params))
}

getClustFreeExprShift <- function(cao.obj, n.top.genes=1000, params=NULL, ...){
  res <- cao.obj$estimateClusterFreeExpressionShifts(n.top.genes=n.top.genes, ...)
  dist.per.type <- res %>% split(cao.obj$sample.groups[names(.)])
  if (is.null(params)) return(dist.per.type)
  return(prepareExpressionDistDf(dist.per.type, params=params))
}

plotExpressionShiftMetric <- function(cao, dist.func, adjust.pvalues=FALSE, n.cores=1, name="expression.shifts", plot.theme=cao$plot.theme,
                                      n.permutations=500, params=NULL, return.all=FALSE, ...) {
  pval.res <- cao$test.results[[name]]$p.dist.info %>%
    estimatePValuesForDistance(sample.groups=cao$sample.groups, diff.func=dist.func, n.cores=n.cores, verbose=FALSE, n.permutations=n.permutations)
  if (adjust.pvalues) pval.res$pvalues %<>% p.adjust(method="BH")
  pval.res$dist.df <- prepareExpressionDistDf(pval.res$dists, params=params)

  if (return.all == "data")
    return(pval.res)

  pval.res$gg <- cacoa:::plotMeanMedValuesPerCellType(pval.res$dist.df, pvalues=pval.res$pvalues, plot.theme=plot.theme, ...)

  if (return.all)
    return(pval.res)

  return(pval.res$gg)
}

generateConsForClustFree <- function(sim, k=30, k.self=5, k.self.weight=0.5, n.cores=30, verbose=TRUE) {
  cms.per.type.per.samp <- sim$cell.groups[colnames(sim$cm)] %>%
    {split(names(.), .)} %>%
    sccore::plapply(function(cg.ids) {
      cm.cg <- sim$cm[,cg.ids]
      sim$sample.per.cell[colnames(cm.cg)] %>% droplevels() %>% {split(names(.), .)} %>%
        lapply(function(sg.ids) cm.cg[,sg.ids])
    }, n.cores=1, progress=verbose)

  p2s.per.type <- cms.per.type.per.samp %>%
    sccore::plapply(lapply, vpscutils::GetPagoda, build.graph=FALSE, verbose=FALSE, n.cores=1, progress=verbose)

  if ((length(p2s.per.type) == 2) & ("value" %in% names(p2s.per.type))) {
    p2s.per.type <- p2s.per.type$value
  }

  con.per.type <- lapply(p2s.per.type, function(p2s) {
    con <- conos::Conos$new(p2s, n.cores=n.cores)
    con$buildGraph(k=k, k.self=k.self, k.self.weight=k.self.weight)
    con
  })

  return(list(p2s.per.type=p2s.per.type, con.per.type=con.per.type, sim=sim))
}

### Sensitivity estimation

prepareSensitivityDf <- function(metric.res, params, covar.name="lfc", trim=0.1) {
  df <- metric.res$dist.df %>% group_by(Type) %>% summarise(value=mean(value, trim=trim)) %>%
    inner_join(params, by=c("Type"="cluster.id")) %>%
    mutate(lfc=factor(paste(lfc), levels=paste(sort(unique(lfc)))),
           ns=factor(paste(ns), levels=paste(sort(unique(ns)))),
           nc=factor(paste(nc), levels=paste(sort(unique(nc)))))

  df$covar <- df[[covar.name]]

  if (!is.null(metric.res$pvalues)) df %<>% mutate(pvalue=metric.res$pvalues[Type])
  return(df)
}

plotSensitivityComparison <- function(sens.df, y=c("distance", "pvalue"), yline=NULL, color.title="LFC",
                                      plot.theme=theme_get(), alpha=0.75, size=2, jitter.width=0.005) {
  y <- match.arg(y)
  if (y == "distance") {
    p.aes <- aes(x=de.frac, y=value, color=covar, shape=ns)
    y.lab <- "Expression distance"
    leg.pos <- c(0, 1)
  } else {
    p.aes <- aes(x=de.frac, y=-log10(pvalue), color=covar, shape=ns)
    y.lab <- "-log10(p-value)"
    leg.pos <- c(1, 0)
    yline <- -log10(0.05)
  }
  gg <- ggplot(sens.df, p.aes) +
    geom_smooth(aes(linetype=ns), size=size, alpha=alpha, se=FALSE) +
    ggbeeswarm::geom_quasirandom(aes(fill=as.factor(rep)), size=size*1.5, alpha=alpha, width=jitter.width) +
    plot.theme + cacoa:::theme_legend_position(leg.pos) +
    labs(x="DE fraction", y=y.lab) +
    theme(legend.background=element_blank()) +
    guides(color=guide_legend(title=color.title), fill=guide_none())

  if (length(unique(sens.df$ns)) < 2) {
    gg <- gg + guides(shape=guide_none(), linetype=guide_none())
  }

  if (!is.null(yline)) gg <- gg + geom_hline(yintercept=yline)
  return(gg)
}

plotSensitivityComparisonPanel <- function(sens.df, yline=NULL, adj.list=NULL, ...) {
  cowplot::plot_grid(
    plotSensitivityComparison(sens.df, y="distance", yline=yline, ...) + adj.list,
    plotSensitivityComparison(sens.df, y="pvalue", ...) + adj.list,
    ncol=2
  )
}

plotExpressionShiftResults <- function(dists.per.type, pvalues, params, ns.col=NULL, trim=0.1, covar.name="lfc", ...) {
  p.df <- list(dist.df=prepareExpressionDistDf(dists.per.type, params=params), pvalues=pvalues) %>%
    prepareSensitivityDf(params=params, covar.name=covar.name, trim=trim)

  if (!is.null(ns.col)) {
    p.df$ns <- p.df[[ns.col]]
  }

  if (is.null(pvalues))
    return(plotSensitivityComparison(p.df, y="distance", yline=0.0, color.title=covar.name, ...))

  return(plotSensitivityComparisonPanel(p.df, yline=0.0, color.title=covar.name, ...))
}

extractRawDistances <- function(p.dist.info, sample.groups) {
  lapply(p.dist.info, function(dm) {
    sample.groups %>% {outer(.[rownames(dm)], .[colnames(dm)], '!=')} %>% dm[.]
  })
}

### Simulating pairwise distances

generateDistanceMatrix <- function(n.cnt=5, n.disease=5, mean.offset=0,
                                   std.cnt=1, std.disease=1, dims=20, do.emb=FALSE) { # change.frac=0.2, dims=500
  cnt.expr <- sapply(1:dims, function(i) rnorm(n.cnt, sd=std.cnt))
  disease.expr <- sapply(1:dims, function(i) rnorm(n.disease, mean=mean.offset, sd=std.disease))
  # disease.expr[,1:round(ncol(disease.expr) * change.frac)] %<>% {. + mean.offset}
  sample.groups <- setNames(rep("A", n.cnt), paste0("A", 1:n.cnt)) %>%
    c(setNames(rep("B", n.disease), paste0("B", 1:n.disease))) %>%
    as.factor()

  expr <- rbind(cnt.expr, disease.expr) %>%
    set_rownames(names(sample.groups)) %>%
    set_colnames(paste0("g", 1:ncol(.)))

  # p.dists <- 1 - cor(t(expr))
  p.dists <- as.matrix(dist(expr))
  emb <- if (do.emb) cmdscale(p.dists, eig=TRUE, k=2)$points else NULL

  return(list(expr=expr, sample.groups=sample.groups, emb=emb, cnt.expr=cnt.expr, disease.expr=disease.expr, p.dists=p.dists))
}

generateDistanceMatrixSimulations <- function(n.samples, mean.offset=0, std.cnt=1, std.disease=1, n.repeats=50) {
  type.labels <- c(both.cross="Shift", ref.cross="Total", ref.target="Variance")
  lapply(names(type.labels), function(n.type) {
    nt <- strsplit(n.type, "\\.")[[1]][1]
    rt <- strsplit(n.type, "\\.")[[1]][2]
    dists <- sapply(1:n.repeats, function(rep) {
      generateDistanceMatrix(n.cnt=n.samples, n.disease=n.samples, mean.offset=mean.offset,
                             std.cnt=std.cnt, std.disease=std.disease) %$%
        cacoa:::estimateExpressionShiftsByDistMat(
          p.dists, sample.groups=sample.groups, ref.level="A", norm.type=nt, return.type=rt
        ) %>% median()
    })

    data.frame(mean.offset=mean.offset, std.cnt=std.cnt, std.disease=std.disease,
               dist=dists, norm.type=type.labels[n.type], n.cnt=n.samples, n.disease=n.samples)
  }) %>% do.call(rbind, .)
}

generateConsForClustFree <- function(sim, k=30, k.self=5, k.self.weight=0.5, n.cores=30, verbose=TRUE, ncomps=30,
                                     metric="angular") {
  cms.per.type.per.samp <- sim$cell.groups[colnames(sim$cm)] %>%
    {split(names(.), .)} %>%
    sccore::plapply(function(cg.ids) {
      cm.cg <- sim$cm[,cg.ids]
      sim$sample.per.cell[colnames(cm.cg)] %>% droplevels() %>% {split(names(.), .)} %>%
        lapply(function(sg.ids) cm.cg[,sg.ids])
    }, n.cores=1, progress=verbose)

  p2s.per.type <- cms.per.type.per.samp %>%
    sccore::plapply(lapply, createPagoda, progress=verbose, n.cores=n.cores)

  if ((length(p2s.per.type) == 2) & ("value" %in% names(p2s.per.type))) {
    p2s.per.type <- p2s.per.type$value
  }

  con.per.type <- sccore::plapply(p2s.per.type, function(p2s) {
    con <- conos::Conos$new(p2s, n.cores=1, verbose=FALSE)
    con$buildGraph(k=k, k.self=k.self, k.self.weight=k.self.weight, verbose=FALSE, ncomps=ncomps, metric=metric)
    con
  }, progress=verbose, n.cores=n.cores)

  return(list(p2s.per.type=p2s.per.type, con.per.type=con.per.type, sim=sim))
}

generateCacoaFromConsForClustFree <- function(con.per.type, sim, n.cores=1) {
  sccore::plapply(con.per.type, function(con) {
    c.ids <- names(con$getDatasetPerCell())
    co <- cacoa::Cacoa$new(
      con, sample.groups=sim$sample.groups[names(con$samples)],
      cell.groups=sim$cell.groups[c.ids], sample.per.cell=droplevels(sim$sample.per.cell[c.ids]),
      ref.level='A', target.level='B', n.cores=n.cores, embedding=NULL
    )

    co$estimateClusterFreeExpressionShifts(gene.selection="expression", verbose=FALSE)
    co
  }, n.cores=1, progress=TRUE)
}

plotClusterFreeShiftSimulations <- function(cao.per.type, params, x.col, adj.list=list()) {
  mapply(function(n, yl) {
    p.df <- cao.per.type %>%
      sapply(function(cao) median(cao$test.results$cluster.free.expr.shifts[[n]], na.rm=TRUE)) %>%
      prepareExpressionDistDf(params=params) %>%
      inner_join(params, c("Type"="cluster.id")) %>%
      mutate(ns=factor(paste0(ns), levels=paste0(sort(unique(ns)))))
    p.df[[x.col]] %<>% factor(levels=sort(unique(.)))

    yline <- if (n == "z_scores") qnorm(0.95) else 0
    ggplot(p.df, aes_string(x=x.col, y='value')) +
      geom_boxplot(outlier.alpha=0, notch=TRUE) +
      ggbeeswarm::geom_quasirandom(size=0.25) +
      geom_hline(yintercept=yline) +
      ylab(yl) +
      cacoa:::theme_legend_position(c(0, 1)) + theme(legend.background=element_blank()) +
      guides(fill=guide_none()) +
      adj.list
  }, c("shifts", "z_scores"), c("Expression distance", "Z-score"), SIMPLIFY=FALSE) %>%
    cowplot::plot_grid(plotlist=., nrow=1)
}

## Boxplot functions

prepareExpressionShiftSimDf <- function(es.res, params=sims$params, sample.groups=sims$sample.groups, sims=NULL) {
  raw.dists <- es.res$p.dist.info %>% extractRawDistances(sample.groups) %>% sapply(median)
  df <- es.res$dists.per.type %>% sapply(median) %>%
    {tibble(Type=names(.), NormDist=., RawDist=raw.dists[names(.)], pvalue=es.res$pvalues[names(.)])} %>%
    inner_join(params, by=c("Type"="cluster.id"))
  return(df)
}

plotExpressionShiftSimDf <- function(p.df, x.col, norm.dist=FALSE, show.pvalues=TRUE, plot.theme=theme_get(),
                                     covar.col=NULL, covar.title=covar.col, dist.col=NULL, adj.list=NULL,
                                     build.panel=TRUE, size=0.1, dodge.width=0.75, alpha=1) {
  if (is.null(p.df$pvalue) || !norm.dist) show.pvalues <- FALSE
  if (is.null(dist.col)) {
    if (norm.dist) {
      dist.col <- 'NormDist'
      y.lab <- 'Normalized distance'
    } else {
      dist.col <- 'RawDist'
      y.lab <- 'Raw distance'
    }
  } else {
    y.lab <- dist.col
  }

  if (!is.null(covar.col)) {
    aes.dist <- aes_string(x=x.col, y=dist.col, fill=covar.col)
    aes.pval <- aes_string(x=x.col, y='pvalue', fill=covar.col)
  } else {
    aes.dist <- aes_string(x=x.col, y=dist.col)
    aes.pval <- aes_string(x=x.col, y='pvalue')
  }


  geom.qr <- ggbeeswarm::geom_quasirandom(aes_string(group=covar.col), size=size, dodge.width=dodge.width, alpha=alpha)
  lg.guide <- guides(fill=guide_legend(title=covar.title))
  fill.scale <- scale_fill_brewer(palette="Dark2")

  gg.dist <- ggplot(p.df, aes.dist) +
    geom_boxplot(notch=TRUE, outlier.alpha=0) +
    geom.qr +
    geom_hline(yintercept=0) +
    scale_y_continuous(expand=c(0, 0, 0.05, 0)) +
    lg.guide + fill.scale + ylab(y.lab) + adj.list

  if (!show.pvalues)
    return(gg.dist)

  gg.pval <- ggplot(p.df, aes.pval) +
    geom_boxplot(notch=TRUE, outlier.alpha=0) +
    geom.qr +
    scale_y_log10() +
    geom_hline(yintercept=0.05) +
    lg.guide + fill.scale + ylab("P-value") + adj.list

  if (build.panel)
    return(cowplot::plot_grid(gg.dist, gg.pval, ncol=2))
  return(list(dist=gg.dist, pval=gg.pval))
}
