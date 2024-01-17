# Mixscale 
Mixscale is an R package designed to analyze CRISPR interference (CRISPRi) based Perturb-seq data. It can quantify the heterogeneity of perturbation strength in each cell and improve the statistical power when doing differential expression (DE) analysis. It also provides functions for downstream analyses including decomposition, permutation test, gene set enrichment test, etc. A brief vignette is available at https://longmanz.github.io/Mixscale/ .

## Dependencies
This package depends on several other R packages:
```
install.packages("Seurat")
install.packages("PMA")
install.packages("protoclust")
BiocManager::install("glmGamPoi")
```

## Installation 
You can easily install the package by the following command:
```
install_github("longmanz/Mixscale")
```
