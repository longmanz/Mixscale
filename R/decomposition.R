#' PCA, MinMax clustering, and MultiCCA analyses, as well as the PCA-based permutation test.


#' The PCA-based permutation test. This function will load in a DE test Z-score
#' matrix, and perform PCA to it to get the 1st to k-th PCs. Then it will permuate
#' the matrix and then the same PCA analysis will be performed. By default, this process
#' will be repeated for 500 times, so that we will have a decent number of null values 
#' in the top k PCs. The proportion of extreme values greater or smaller than the actual
#' value will be used as the P-values. 
#' 
#' @export
#' 
#' @param mat the Z-score matrix to perform the permutation test. 
#' Rows are the gene and columns are the conditions/samples.
#' @param k the number of top PCs to extract. If set to NULL, the 
#' function will determine the number of PCs to extract based on 
#' the var_prop, which is the cut-off on the proportion of variance 
#' explained for each PC. Only PCs with prop_var greater than this 
#' value will be removed extracted.
#' @param var_prop if k is not set, then use this value as a cutoff for 
#' %var explained to select the PCs. Only the PCs with %var larger than this 
#' value will be selected (the top k PCs). 
#' @param center a boolen value to indicate whether the column will be centered.
#' @param scale a boolen value to indicate whether the column will be scaled to 1.
#' @param num_iter the number of iteration for the permutation test.
#' 
#' @param seed seed for random number generator.
#' @return a list object consists of 3 elements: the original input Z-score matrix,
#' the p-value matrix produced by the permutation test for each genes (same dimension 
#' as the original input matrix), and the prcomp() object done for the input Z-score 
#' matrix.
#' 

PCApermtest = function(mat = NULL, k = 1, var_prop = 0.1, 
                       center = T, scale = T, 
                       num_iter = 200, 
                       seed = 123124125,
                       ...){
    # set the random seed 
    set.seed(seed = seed)
    
    # perform standard PCA (centering and scaling) for the input mat
    test = prcomp(mat, center = center, scale. = scale)
    
    # get the % of variance (which is sdev^2) explained for each PC 
    prop_test = test$sdev^2/sum(test$sdev^2)
    
    # if k is not set, we will determine it based on the %var by each PC
    if(is.null(k) & !is.null(var_prop)){
        # we will only use those PCs with %var >= 10% 
        k = tail(which(prop_test >= var_prop), 1)
        # check if none of the PC has %var >= prop_var.
        if(length(k) == 0){
            stop(paste("None of the PC has %var >=", var_prop, ". Please use a lower value."))
        }
    } else if (is.null(k) & is.null(var_prop)){
        stop(paste("Neither k or var_prop is provided. Please specify one of them."))
    }
    
    # the permutation test 
    # Preallocate memory for null_pc
    null_pc <- matrix(0, nrow = nrow(mat) * num_iter, ncol = k)
    
    # Perform iterations
    for(idx_iter in 1:num_iter){
        null_tmp = apply(mat, 2, sample)
        null_test = prcomp(null_tmp, center = center, scale. = scale)
        
        # Direct indexing to assign values
        start_idx <- (idx_iter - 1) * nrow(mat) + 1
        end_idx <- idx_iter * nrow(mat)
        
        null_pc[start_idx:end_idx, ] <- null_test$x[1:k]
    }
    
    #### after getting the null distribution we will now move on to pval calculation
    for(i in 1:ncol(null_pc)){
        ecdf_fun = ecdf(null_pc[, i])
        tmp_pval = sapply(X = test$x[, i], FUN = ecdf_fun)
        if(i == 1){
            pval = tmp_pval
        } else {
            pval = cbind(pval, tmp_pval)
        }
    }
    colnames(pval) = paste0("PC", 1:ncol(null_pc))
    
    return(list(mat = mat, 
                pmat = pval, 
                prcomp_obj = test))

}


#' The function will load the pval matrix calculated from PCApermtest() and 
#' return the significant rows (usually genes) given the threshold that users
#' provide. 
#' 
#' @export
#' @param perm_obj the list object produces by the PCApermtest() function.
#' @param k the number of top PCs to extract. If set to NULL, the 
#' function will determine the number of PCs to extract based on 
#' the var_prop, which is the cut-off on the proportion of variance 
#' explained for each PC. Only PCs with prop_var greater than this 
#' value will be removed extracted.
#' @param var_prop if k is not set, then use this value as a cutoff for 
#' %var explained to select the PCs. Only the PCs with %var larger than this 
#' value will be selected (the top k PCs). 
#' @param var_prop_total similar to var_prop. The accumulated sum of the %var
#' of the top i = 1, 2, ... PCs will be calculated, and the top k PCs with the 
#' accumulated sum larger than this value to be selected. 
#' @param perm_pval_thres the p-value threshold for the permutation test. Rows 
#' (genes) with permutation p-value lower then permtest_pval_thres or greater than 
#' (1 - permtest_pval_thres) will be selected. Default is 0.05.
#' @param ori_pval_thres the p-value threshold for the original DE test. In the original
#' Z-score matrix we store all the Z-score from the DE test. When selecting significant 
#' rows (genes), we might wish them to be significant in the original DE test as well. 
#' By setting this value, we force the rows (genes) to be both significant in the 
#' permutation test and the DE test. Default is 0.05/30000 = 1.666667e-06 (which is 
#' gene-wide significant after Bonferroni correction). Set it to 1 to avoid such filtering.
#' @param cor_threshold After the PCApermtest, the actual orientation of the PC
#' might not be the same as the orientation of its correlated columns in the original
#' matrix. We need to do correlation test between each PC and all the columns in the 
#' original matrix, and then use this cor_threshold to define and extract the 
#' correlated columns and then use them to determine the actual orientation of the PC.
#' @param collapse a boolen value to indicate when k >= 2, whether the significant 
#' genes from the top k PCs should be return together as one gene list (collapse = T) 
#' or separately for each k (collapse = F). Default is T.
#' 
#' @return a list of significant genes (gene IDs) selected from the output of PCApermtest. 
#'  The order of genes in each gene list indicates the significance of P-value (high
#'  significance at the top).

get_sig_genes = function(perm_obj = NULL, 
                         k = 1, 
                         var_prop = NULL,
                         var_prop_total = NULL,
                         perm_pval_thres = 0.05, 
                         ori_pval_thres = 1.666667e-06, 
                         cor_threshold = 0.2, 
                         collapse = T, 
                         ...){
    # check if the object has the correct format
    if(class(perm_obj$prcomp_obj) != "prcomp"){
        stop("The input object does not contain the correct PCA obj from prcomp().")
    }
    
    # extract each element from the provided object 
    test = perm_obj$prcomp_obj
    mat = perm_obj$mat
    pval = perm_obj$pmat
    
    # first we need to re-calculate the %var explained for each PC
    prop_test = test$sdev^2/sum(test$sdev^2)
    
    # if k is set then we will use the top k PCs.
    # if not, we will determine the k based on the given var_prop or var_prop_total
    #  note that if var_prop is set, var_prop_total will be ignored. 
    if(!is.null(k)){
        max_PC_idx = k
    } else if (is.null(k) & !is.null(var_prop)) {
        idx1 = which(prop_test >= var_prop)
        max_PC_idx = tail(idx1, 1)
    } else if (is.null(k) & is.null(var_prop) & !is.null(var_prop_total)){
        max_PC_idx = which( cumsum(prop_test) >= var_prop_total )[1]
    } else {
        stop("Please provide at least one of the following paramters: k, var_prop, var_prop_total.")
    }
    
    # do the selection
    slct_pval = pval[, max_PC_idx, drop = F]
    
    # extract the sig elements for each PC
    z_threshold = sqrt(qchisq(ori_pval_thres, df=1, lower=F))
    # 
    for(i in 1:max_PC_idx){
        topDEG_idx = list()
        bottomDEG_idx = list()
        
        cor_test = cor(test$x[,i], mat)
        print(cor_test)
        PRTB_idx = which(abs(cor_test) >= cor_threshold)
        
        topDEG = which(slct_pval[, i] >= 1 - pval_threshold & apply(mat[, PRTB_idx, drop=F], MARGIN = 1, FUN = function(x) any(abs(x) >= z_threshold) ))
        topDEG = topDEG[order(test$x[topDEG, i], decreasing = T)]
        bottomDEG = which(slct_pval[, i] <= pval_threshold & apply(mat[, PRTB_idx, drop=F], MARGIN = 1, FUN = function(x) any(abs(x) >= z_threshold) ))
        bottomDEG = bottomDEG[order(test$x[bottomDEG, i], decreasing = T)]
        
        # determine the orientation by looking at the sums of Z-score
        # a positive sum(Z-score) indicates up-regulated genes; 
        # a negative sum(Z-score) indicates down-regulated genes.
        if( sum(mat[topDEG, PRTB_idx]) > sum(mat[bottomDEG, PRTB_idx]) ){
            swap = topDEG
            topDEG = rev(bottomDEG)
            bottomDEG = rev(swap)
            rm(swap)
        }
        topDEG_idx[[paste0("PC", i)]] = topDEG
        bottomDEG_idx[[paste0("PC", i)]] = bottomDEG
        
    }
    
    # By default, if k >= 2, we will collapse all the significant genes together
    #  and return it to users. 
    
    if(collapse == T){
        topDEG_idx = Reduce(c, topDEG_idx)
        topDEG_idx = topDEG_idx[!duplicated(topDEG_idx)]
        
        bottomDEG_idx = rev(bottomDEG_idx)
        bottomDEG_idx = Reduce(c, bottomDEG_idx)
        bottomDEG_idx = bottomDEG_idx[!duplicated(bottomDEG_idx)]
        
        return(list(downDEGs = rownames(mat)[topDEG_idx], 
                    upDEGs = rownames(mat)[bottomDEG_idx]))
    } else {
        downDEGs = list()
        upDEGs = list()
        for(i in 1:max_PC_idx){
            downDEGs[[names(topDEG_idx)[i]]] = rownames(mat)[topDEG_idx[[1]]]
            upDEGs[[names(bottomDEG_idx)[i]]] = rev(rownames(mat)[bottomDEG_idx[[1]]])
        }
        return(list(downDEGs = downDEGs, 
                    upDEGs = upDEGs))
    }

}



#' A wrapper for different hierarchical clustering methods to be applied to the within-cell-type
#' cross-conditions Z-score matrix. Highly similar conditions (columns) will be grouped together 
#' given the DE Z-scores of rows (genes). 
#' 
#' @param mat the Z-score matrix to perform the permutation test. 
#' Rows are the gene and columns are the conditions/samples.
#' @param cor_method the method to calculate the correlation matrix. see cor() 
#' function for details. 
#' @param hclust_method the method to perform the hierarchical clustering. The default method
#' is MinMax clustering (package 'protoclust' required). Other methods available in the hclust()
#' function is also allowed. 
#' @param dist_thres The distance to cut the hierarchical clustering tree ("tree height"). See cutree() 
#' function for details. Default is 0.6.
#' 
#' @return return a list of two object: 1. a list of the cluster assignment of the 
#' columns (only those got successfully assigned to a multi-member cluster will be stored).
#' 2. a object of the output object from protocut() (if MinMax hclust method is selected) or 
#' from cutree() (if other hclust method is selected)
#' 
#' 

HclustTest = function(mat = NULL, 
                      cor_method = c("pearson", "kendall", "spearman"), 
                      hclust_method = "minmax", 
                      dist_thres = 0.6,
                      ...){

    # skip if the number of column is <= 2
    if(ncol(mat) <= 2){
        stop("The number of columns in the input matrix is less than 3. Hierarchical clustering cannot be applied.")
    }
    
    # calculate the distance matrix ( = 1 - correlation_matrix)
    cor_mat = abs(cor(mat, method = cor_method))
    mat.dist = as.dist(1 - cor_mat)
    
    # run the hclust using the input method
    if(hclust_method == "minmax"){
        # check if the required package is installed or not
        protoclust.installed <- PackageCheck("protoclust", error = FALSE)
        if (!protoclust.installed[1]) {
            stop("Please install the protoclust package to use MinMax hierarchical clustering method.", 
                 "\nThis can be accomplished with the following command: ", 
                 "\n----------------------------------------", "\ninstall.packages('protoclust')", 
                 "\n----------------------------------------", call. = FALSE)
        }
        # 
        mat.tree = protoclust(mat.dist)
        mat.cut <- protocut(mat.tree, h = dist_thres)
        cl = mat.cut$cl  # the cluster assignment
    } else {
        mat.tree <- hclust(mat.dist, method="single")
        mat.cut <- cutree(mat.tree, h = dist_thres)
        cl = mat.cut
    }
    
    # count of elements/points in each cluster, save the IDs of any cluster that has >= 2 members
    ct_pt = table(cl)
    slct_cluster = names(ct_pt)[which(ct_pt > 1)]
    slct_cluster = as.integer(slct_cluster)
    
    #  save the cluster information (members) for those slct_cluster
    cluster_list = list()
    t = 1
    for(num_cluster in slct_cluster){
        item = names(cl)[cl == num_cluster]
        cluster_list[[t]] = item
        t = t + 1
    }
    
    # return the cluster information
    return(list = c(cluster_assignment = cluster_list, 
                    full_obj = mat.cut))

}









