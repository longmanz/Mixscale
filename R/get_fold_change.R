#'
NULL


#' Calculate log-fold-change given a vector of gene expression and the indices of perturbed cells and non-target cells
#'
#' Function to calculate log-fold-change for pooled CRISPR screen datasets.
#' It is just a simple function to calculate the log-fold-change. Users can customise the min.cells, 
#' minimal expression threshold, pseudo-count (the small value added to the expression level to avoid log(0)), 
#' minimal percentage of cells expression the genes, and the base of the log.  
#'
#' @param gene_exp a vector of the gene expression levels
#' @param idx_P a vector of index for the perturbed cells in the gene_exp
#' @param idx_NT a vector of index for the non-target cells (controls) in the gene_exp
#' @param min.cells the minimal number of cells that expresses the gene; if lower than this value, the 
#' fold-change will be returned as NA. Default is 3. 
#' @param thresh.min the minimal value of expression; any expression value lower than this will be 
#' considered as 0. Default is 0.
#' @param pseudocount.use the small value that will be added to the log-transformation to avoid log(0).
#' For example, if a mean expression value is x, the final log-
#' @param min.pct the minimal proportion of cells in either groups that expresses the gene
#' @param base the base for log()
#' @param norm.method the normalization method for the input gene_exp. Default is 'raw', which means 
#' the original count value without normalization. The other supported values are 'log.norm', "scale.data". 
#' The mean.fxn() will change accordingly. 
#' @return Returns a single value of the log-fold-change of the input gene.
#' @export
#' @concept perturbation_scoring

get_fc = function(gene_exp = NULL, idx_P = NULL, idx_NT = NULL, 
                  min.cells = 3,  
                  thresh.min = 0, 
                  pseudocount.use = 1,  
                  min.pct = 0.1,  
                  base = 2, 
                  norm.method = 'raw'
){
    # the exp_vec should have been std
    
    # flag 1: do minimum cell check
    if (sum(gene_exp[idx_P] > 0) < min.cells &&
        sum(gene_exp[idx_NT] > 0) < min.cells) {
        return(NA)
    } 
    
    # flag 2: do variance check (not 0)
    if (var(gene_exp) == 0) {
        return(NA)
    }
    
    # flag 3: min.pct check
    #   calculate fraction of cell with expression > 0
    pct.1 <- round(
        x = sum(x = gene_exp[idx_NT] > thresh.min) /
            length(x = gene_exp[idx_NT]),
        digits = 3
    )
    pct.2 <- round(
        x = sum(x = gene_exp[idx_P] > thresh.min) /
            length(x = gene_exp[idx_P]),
        digits = 3
    )
    
    if (pct.1 < min.pct & pct.2 < min.pct) {
        return(NA)
    } 
    
    # set the mean.fxn() function according to norm.method
    default.mean.fxn <- function(x) {
        return(log(x = mean(x = x) + pseudocount.use, base = base))
    }
    mean.fxn <- switch(
        EXPR = norm.method,
        'log.norm' = function(x) {
                return(log(x = mean(x = expm1(x = x)) + pseudocount.use, base = base))
            },
        'scale.data' = mean,
        default.mean.fxn
    )
    
    # flag 4: minimum fold change check
    # log(x = mean(x = x) + pseudocount.use, base = base)
    data.1 <- mean.fxn(gene_exp[idx_NT])
    data.2 <- mean.fxn(gene_exp[idx_P])
    fc <- - data.1 + data.2
    return(fc)
    
}



# a function to do the filtering for all the genens
get_idx = function(gene_exp = NULL, idx_P = NULL, idx_NT = NULL, 
                   min.cells = 3,  # the minimum cell threshold to perform DE
                   thresh.min = 0, # the minimum expression level
                   pseudocount.use = 1,  
                   min.pct = 0.1,  
                   logfc.threshold = 0.1,
                   base = 2, 
                   norm.method = 'raw'
){

    fc <- get_fc(gene_exp = gene_exp, idx_P = idx_P, idx_NT = idx_NT, 
                 min.cells = min.cells, 
                 thresh.min = thresh.min,
                 pseudocount.use = pseudocount.use,  
                 min.pct = min.pct,  
                 base = base, 
                 norm.method = norm.method
                 )
    
    if (fc < logfc.threshold){
        return(FALSE)
    }
    
    return(TRUE)
}


#' Calculate log-fold-change given a vector of gene expression and the indices of perturbed cells and non-target cells
#'
#' Function to calculate log-fold-change for pooled CRISPR screen datasets.
#' It is just a simple function to calculate the log-fold-change. Users can customise the min.cells, 
#' minimal expression threshold, pseudo-count (the small value added to the expression level to avoid log(0)), 
#' minimal percentage of cells expression the genes, and the base of the log.  
#'
#' @inheritParams Seurat::FoldChange
#' @return Returns a single value of the log-fold-change of the input gene.
#' @export
#' @concept perturbation_scoring
FoldChange_new <- function(
        object,
        cells.1,
        cells.2,
        mean.fxn,
        fc.name,
        features = NULL,
        ...
) {
    features <- features %||% rownames(x = object)
    
    # Calculate percent expressed
    thresh.min <- 0
    
    min.cell.1 = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) 
    min.cell.2 = rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) 
    
    pct.1 <- round(
        x = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
            length(x = cells.1),
        digits = 3
    )
    pct.2 <- round(
        x = rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
            length(x = cells.2),
        digits = 3
    )
    # Calculate fold change
    data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
    data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
    fc <- (data.1 - data.2)
    fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2, min.cell.1, min.cell.2))
    colnames(fc.results) <- c(fc.name, "pct.1", "pct.2", "min.cell.1", "min.cell.2")
    return(fc.results)
}




