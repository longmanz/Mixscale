# Perturbation_Scoring

## Dependencies
This package depends on several other R packages:
```
install.packages("Seurat")
install.packages("PMA")
install.packages("ggplot2")
install.packages("ggridges")
install.packages("protoclust")
BiocManager::install("glmGamPoi")
```

## Installation 
You can easily install the package by the following command:
```
install_github("longmanz/PRTBScoring")
```
\
Or, if you prefer to download it first and then install it locally, please use the following command
```
git clone https://github.com/longmanz/PRTBScoring.git
devtools::load_all("/path_to_the_downloaded_directory/")
```

