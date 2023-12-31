% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomposition.R
\name{DEhclust}
\alias{DEhclust}
\title{Run Hierarchical clustering for a matrix}
\usage{
DEhclust(
  mat = NULL,
  cor_method = c("pearson", "kendall", "spearman"),
  hclust_method = "minmax",
  dist_thres = 0.6,
  ...
)
}
\arguments{
\item{mat}{the Z-score matrix to perform the permutation test. 
Rows are the gene and columns are the conditions/samples.}

\item{cor_method}{the method to calculate the correlation matrix. see cor() 
function for details.}

\item{hclust_method}{the method to perform the hierarchical clustering. The default method
is MinMax clustering (package 'protoclust' required). Other methods available in the hclust()
function is also allowed.}

\item{dist_thres}{The distance to cut the hierarchical clustering tree ("tree height"). See cutree() 
function for details. Default is 0.6.}
}
\value{
return a list of two object: 1. a list of the cluster assignment of the 
columns (only those got successfully assigned to a multi-member cluster will be stored).
2. a object of the output object from protocut() (if MinMax hclust method is selected) or 
from cutree() (if other hclust method is selected)
}
\description{
A wrapper for different hierarchical clustering methods to be applied to the within-cell-type
cross-conditions Z-score matrix (input). Highly similar conditions (columns) will be grouped together 
given the DE Z-scores of rows (genes).
}
