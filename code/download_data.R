library(dataorganizer)
library(magrittr)

geoUrl <- function(...) paste0('https://www.ncbi.nlm.nih.gov/geo/download/?format=file&', ...)

downloadFileToData <- function(url, path, force=FALSE, gunzip=FALSE, unzip=FALSE) {
  path <- DataPath(path)
  if (!force && file.exists(path))
    return()

  if (gunzip) path <- paste0(path, ".gz")
  download.file(url, path)
  if (gunzip) {
    R.utils::gunzip(path)
  } else if (unzip) {
    utils::unzip(path, exdir=dirname(path))
  }
}

downloadData <- function(datasets=NULL) {
  all.datasets <- c("AZ", "ASD", "EP", "MS", "PF", "SCC")

  if (is.null(datasets)) datasets <- all.datasets
  if (length(setdiff(datasets, all.datasets)) > 0)
    stop("Unknown dataset: ", paste0(setdiff(datasets, all.datasets), collapse=','))

  for (ds in datasets) {
    message(ds)
    suppressWarnings(dir.create(DataPath(ds)))
    if (ds == "AZ") {
      # Metadata was downloaded from http://adsn.ddnetbio.com/ , but it doesn't allow persistent links, so we copied it to our server
      downloadFileToData('http://pklab.med.harvard.edu/viktor/publications/Cacoa/alzheimer_metadata.tsv', "AZ/cell_metadata.csv")
      downloadFileToData(geoUrl('acc=GSE138852&file=GSE138852%5Fcounts%2Ecsv%2Egz'), "AZ/cell_counts.csv", gunzip=TRUE)
    } else if (ds == "ASD") {
      downloadFileToData("https://cells.ucsc.edu/autism/rawMatrix.zip", "ASD/rawMatrix.zip", unzip=TRUE)
    } else if (ds == "EP") {
      downloadFileToData("https://raw.githubusercontent.com/khodosevichlab/Epilepsy19/master/metadata/annotation.csv", "EP/annotation.csv")
      downloadFileToData("https://raw.githubusercontent.com/khodosevichlab/Epilepsy19/master/metadata/sample_info.csv", "EP/sample_info.csv")
      downloadFileToData("http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_samples.rds", "EP/con_filt_samples.rds")
    } else if (ds == "MS") {
      downloadFileToData("http://cells.ucsc.edu/ms/rawMatrix.zip", "MS/rawMatrix.zip", unzip=TRUE)
    } else if (ds == "PF") {
      downloadFileToData(geoUrl('acc=GSE135893&file=GSE135893%5Fbarcodes%2Etsv%2Egz'), "PF/barcodes.tsv", gunzip=TRUE)
      downloadFileToData(geoUrl('acc=GSE135893&file=GSE135893%5Fgenes%2Etsv%2Egz'), "PF/genes.tsv", gunzip=TRUE)
      downloadFileToData(geoUrl('acc=GSE135893&file=GSE135893%5Fmatrix%2Emtx%2Egz'), "PF/matrix.mtx", gunzip=TRUE)
      downloadFileToData(geoUrl('acc=GSE135893&file=GSE135893%5FIPF%5Fmetadata%2Ecsv%2Egz'), "PF/cell_metadata.csv", gunzip=TRUE)
    } else if (ds == "SCC") {
      downloadFileToData(geoUrl('acc=GSE144236&file=GSE144236%5FcSCC%5Fcounts%2Etxt%2Egz'), "SCC/counts.txt", gunzip=TRUE)
      downloadFileToData(geoUrl('acc=GSE144236&file=GSE144236%5Fpatient%5Fmetadata%5Fnew%2Etxt%2Egz'), "SCC/cell_metadata.txt", gunzip=TRUE)
    }
  }
}

downloadData()
