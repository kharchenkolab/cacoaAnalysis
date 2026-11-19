# Prepare the cSCC dataset (Ji et al. 2020) for scc_cacoa.Rmd.
# Reads the Cacoa object of the Cacoa preprint (cacoaAnalysis/data/SCC/cao.rds: a Conos object of pagoda2 samples plus the cell annotation and the
# embedding), refreshes the R6 objects to the installed package versions, and stores the pieces as one list:
#   con (Conos), cell.groups (named factor), sample.per.cell (named factor), embedding (cells x 2).
# usage: Rscript prepare_scc_data.R <out.rds> [<old cacoa rds>]
suppressMessages({library(conos); library(pagoda2); library(magrittr)})
args <- commandArgs(trailingOnly = TRUE)
out.file <- args[1]
old.file <- if (length(args) > 1) args[2] else "/home/viktor_petukhov/projects/cacoaAnalysis/data/SCC/cao.rds"   # (the older object of 2021, 18 samples, was /home/viktor_petukhov/data/cacoa/scc.cacoa.rds)

ocao <- readRDS(old.file)
for (n in names(ocao$data.object$samples)) ocao$data.object$samples[[n]] %<>% pagoda2::Pagoda2$new()
con <- Conos$new(ocao$data.object)
con$graph <- igraph::upgrade_graph(con$graph)   # the stored graph predates the installed igraph

# the source count matrices carry two metadata columns as genes (the patient number and the tissue flag, 1 = tumor);
# they would act as perfect patient and condition markers, so they are dropped
META <- c("Patient", "Tissue: 0=Normal, 1=Tumor")
for (n in names(con$samples)) {
  p2 <- con$samples[[n]]
  p2$misc$rawCounts <- p2$misc$rawCounts[, !colnames(p2$misc$rawCounts) %in% META, drop = FALSE]
  if (!is.null(p2$counts)) p2$counts <- p2$counts[, !colnames(p2$counts) %in% META, drop = FALSE]
}

# the old annotation merges the paper's tumour and normal keratinocyte states (KC_Basal / KC_Cyc / KC_Diff) and calls
# the eccrine cells "Epithelial"; the paper's name is used here
cell.groups <- ocao$cell.groups
levels(cell.groups)[levels(cell.groups) == "Epithelial"] <- "Eccrine"

saveRDS(list(con = con, cell.groups = cell.groups, sample.per.cell = ocao$sample.per.cell, embedding = ocao$embedding), out.file)
