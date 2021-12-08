# Data

The repo uses data from the following publications:

1. [Alzheimer](https://www.nature.com/articles/s41593-019-0539-4/): Grubman, A., Chew, G., Ouyang, J. F., Sun, G., Choo, X. Y., McLean, C., ...Polo, J. M. (2019). A single-cell atlas of entorhinal cortex from individuals with Alzheimer’s disease reveals cell-type-specific gene expression regulation - Nature Neuroscience. Nat. Neurosci., 22, 2087–2097. doi: 10.1038/s41593-019-0539-4
2. [Autism](https://www.science.org/doi/10.1126/science.aav8130): Velmeshev, D., Schirmer, L., Jung, D., Haeussler, M., Perez, Y., Mayer, S., ...Kriegstein, A. R. (2019). Single-cell genomics identifies cell type–specific molecular changes in autism. Science.
3. [Epilepsy](https://www.nature.com/articles/s41467-020-18752-7): Pfisterer, U., Petukhov, V., Demharter, S., Meichsner, J., Thompson, J. J., Batiuk, M. Y., ...Khodosevich, K. (2020). Identification of epilepsy-associated neuronal subtypes and gene expression underlying epileptogenesis - Nature Communications. Nat. Commun., 11(5038), 1–19. doi: 10.1038/s41467-020-18752-7
4. [MS](https://www.nature.com/articles/s41586-019-1404-z): Schrimer et al. “Neuronal vulnerability and multilineage diversity in multiple sclerosis”, Nature 2019 Sep;573(7772):75-82. doi: 10.1038/s41586-019-1404-z. ACCESSION ID
5. [PF](https://www.science.org/doi/10.1126/sciadv.aba1972): Habermann, A. C., Gutierrez, A. J., Bui, L. T., Yahn, S. L., Winters, N. I., Calvi, C. L., ...Kropski, J. A. (2020). Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. Sci. Adv. Retrieved from https://www.science.org/doi/10.1126/sciadv.aba1972
6. [SCC](https://www.sciencedirect.com/science/article/pii/S0092867420306723): Ji, A. L., Rubin, A. J., Thrane, K., Jiang, S., Reynolds, D. L., Meyers, R. M., ...Khavari, P. A. (2020). Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma. Cell, 182(2), 497–514.e22. doi: 10.1016/j.cell.2020.05.039

Most notebooks are on Bayes: /home/meisl/Workplace/cacoa/data/Script/
MS notebook on mendel:/d0/home/pkharchenko/p2/comp/MS/examplle.Rmd

## Downloading

To download all data automatically, run `Rscript code/download_data.R` from the project root.
The descriptions below are for downloading by hand and navigating data repositories.

### Alzheimer

[Interactive exploration](http://adsn.ddnetbio.com/)

```
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138852&format=file&file=GSE138852%5Fcovariates%2Ecsv%2Egz -O GSE138852_covariates.csv.gz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138852&format=file&file=GSE138852%5Fcounts%2Ecsv%2Egz -O GSE138852_counts.csv.gz
gunzip ./GSE138852*
```

### Autism

[Interactive exploration](https://cells.ucsc.edu/?ds=autism)

```
wget https://cells.ucsc.edu/autism/rawMatrix.zip
unzip rawMatrix.zip
```

Returns 10x in .mm format and cell metadata in meta.txt

### Epilepsy

[Interactive exploration](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_samples.bin)

```
wget https://raw.githubusercontent.com/khodosevichlab/Epilepsy19/master/metadata/annotation.csv
wget https://raw.githubusercontent.com/khodosevichlab/Epilepsy19/master/metadata/sample_info.csv
wget pklab.med.harvard.edu/viktor/publications/Epilepsy19/count_matrices.rds
```

Optionally:

```
wget http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_samples.rds
```

### Multiple sclerosis

[Ineteractive exploration](http://cells.ucsc.edu/?ds=ms)

```
wget http://cells.ucsc.edu/ms/rawMatrix.zip
```

Same as Autism above.

### Pulmonary fibrosis

```
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135893&format=file&file=GSE135893%5Fbarcodes%2Etsv%2Egz -O GSE135893_barcodes.tsv.gz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135893&format=file&file=GSE135893%5Fgenes%2Etsv%2Egz -O GSE135893_genes.tsv.gz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135893&format=file&file=GSE135893%5Fmatrix%2Emtx%2Egz -O GSE135893_matrix.mtx.gz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135893&format=file&file=GSE135893%5FIPF%5Fmetadata%2Ecsv%2Egz -O GSE135893_IPF_metadata.csv.gz
gunzip GSE135893_*
```

### Squamous Cell Carcinoma

If it's a wrong file, need to check https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144236

```
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144236&format=file&file=GSE144236%5FcSCC%5Fcounts%2Etxt%2Egz -O GSE144236_cSCC_counts.txt.gz
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144236&format=file&file=GSE144236%5Fpatient%5Fmetadata%5Fnew%2Etxt%2Egz -O GSE144236_patient_metadata_new.txt.gz
```
