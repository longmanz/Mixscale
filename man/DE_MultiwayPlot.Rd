% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{DE_MultiwayPlot}
\alias{DE_MultiwayPlot}
\title{Multi-way dotplot for multi-cell-line DE results}
\usage{
DE_MultiwayPlot(DEG_mat = NULL, zscore_cap = 10, ...)
}
\arguments{
\item{name}{description}
}
\description{
This function will generate a multi-way dotplot if the DE results produces by scoringDE() 
and get_DE_mat() contain multiple different cell-lines. The a-axis will be some selected perturbations and 
the y-axis will be some selected DE genes that users want to explore. Within each column, 
multiple dots will be displayed, with size representing the DE test DE Z-score (significance)
and color indicating the cell line identity. It is a good way to explore the consistency and
heterogeneity of DE results across perturbations/cell lines.
}
