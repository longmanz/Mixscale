#'
NULL

#' Run PCA-based permutation test for a matrix
#' 
#' This function will load in a DE test Z-score matrix, and perform PCA to it to get the 1st to k-th PCs. 
#' Then it will permuate the matrix and then the same PCA analysis will be performed. By default, this process
#' will be repeated for 200 times (default), so that we will have a decent number of null values 
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
#' @param var_prop_total similar to var_prop. The accumulated sum of the %var
#' of the top i = 1, 2, ... PCs will be calculated, and the top k PCs with the 
#' accumulated sum larger than this value to be selected. 
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

PCApermtest = function(mat = NULL, k = 1, 
                       var_prop = NULL, 
                       var_prop_total = NULL,
                       center = T, scale = T, 
                       row_filtering_pval = 0.05,
                       num_iter = 200, 
                       seed = 123124125,
                       ...){
    # set the random seed 
    set.seed(seed = seed)
    
    # filter the row if row_filtering_pval is set:
    if(!is.null(row_filtering_pval) ){
        if(row_filtering_pval > 0 & row_filtering_pval <= 1){
            rm_row_idx = which( apply(mat, MARGIN = 1, FUN = function(x) all(abs(x) < sqrt(qchisq(row_filtering_pval, df = 1, lower.tail = F)) )) )
            print(paste("Removing", length(rm_row_idx), "rows given the row_filtering_pval =", row_filtering_pval))
            # remove rows if there is any un-qualified row
            if(length(rm_row_idx) != 0){
                mat = mat[-rm_row_idx, , drop = F]
            }
        } else {
            print("Please set a valid row_filtering_pval that is between 0 and 1.")
        }
    } 
    
    # perform standard PCA (centering and scaling) for the input mat
    test = prcomp(mat, center = center, scale. = scale)
    
    # get the % of variance (which is sdev^2) explained for each PC 
    prop_test = test$sdev^2/sum(test$sdev^2)
    
    # if k is not set, we will determine it based on the %var by each PC
    if(!is.null(k)){
        # 
    } else if (is.null(k) & !is.null(var_prop)) {
        k = tail(which(prop_test >= var_prop), 1)
        if(length(k) == 0){
            stop(paste("None of the PC has %var >=", var_prop, ". Please use a lower value."))
        }
    } else if (is.null(k) & is.null(var_prop) & !is.null(var_prop_total)){
        k = which( cumsum(prop_test) >= var_prop_total )[1]
    } else {
        stop("Please provide at least one of the following paramters: k, var_prop, var_prop_total.")
    }
    
    
    # the permutation test 
    # Preallocate memory for null_pc
    null_pc <- matrix(0, nrow = nrow(mat) * num_iter, ncol = k)
    
    # Perform iterations
    for(idx_iter in 1:num_iter){
        # set.seed(idx_iter+500)
        
        null_tmp = apply(mat, 2, sample)
        null_test = prcomp(null_tmp, center = center, scale. = scale)
        
        # Direct indexing to assign values
        start_idx <- (idx_iter - 1) * nrow(mat) + 1
        end_idx <- idx_iter * nrow(mat)
        
        null_pc[start_idx:end_idx, ] <- null_test$x[, 1:k, drop = F]
    }
    
    #### after getting the null distribution we will now move on to pval calculation
    for(i in 1:ncol(null_pc)){
        ecdf_fun = ecdf(null_pc[, i])
        tmp_pval = sapply(X = test$x[, i], FUN = ecdf_fun)
        if(i == 1){
            pval = matrix(data = tmp_pval, ncol = 1)
        } else {
            pval = cbind(pval, tmp_pval)
        }
    }
    colnames(pval) = paste0("PC", 1:ncol(null_pc))
    rownames(pval) = rownames(mat)
    
    return(list(mat = mat, 
                pmat = pval, 
                prcomp_obj = test))

}


#' Extract significant genes from PCApermtest
#' 
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
    slct_pval = pval[, 1:max_PC_idx, drop = F]
    
    # extract the sig elements for each PC
    # note that pval_threshold is not the same as perm_pval_thres
    pval_threshold = ori_pval_thres
    z_threshold = sqrt(qchisq(pval_threshold, df=1, lower=F))
    
    ### 
    topDEG_idx = list()
    bottomDEG_idx = list()
    for(i in 1:max_PC_idx){

        cor_test = cor(test$x[,i], mat)
        # print(cor_test)
        PRTB_idx = which(abs(cor_test) >= cor_threshold)
        
        topDEG = which(slct_pval[, i] >= 1 - perm_pval_thres & apply(mat[, PRTB_idx, drop=F], MARGIN = 1, FUN = function(x) any(abs(x) >= z_threshold) ))
        topDEG = topDEG[order(test$x[topDEG, i], decreasing = T)]
        bottomDEG = which(slct_pval[, i] <= perm_pval_thres & apply(mat[, PRTB_idx, drop=F], MARGIN = 1, FUN = function(x) any(abs(x) >= z_threshold) ))
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
        
        bottomDEG_idx = lapply(bottomDEG_idx, rev)
        bottomDEG_idx = Reduce(c, bottomDEG_idx)
        # bottomDEG_idx = rev(bottomDEG_idx)
        bottomDEG_idx = bottomDEG_idx[!duplicated(bottomDEG_idx)]
        
        return(list(downDEGs = rownames(mat)[topDEG_idx], 
                    upDEGs = rownames(mat)[bottomDEG_idx]))
    } else {
        downDEGs = list()
        upDEGs = list()
        for(i in 1:max_PC_idx){
            downDEGs[[names(topDEG_idx)[i]]] = rownames(mat)[topDEG_idx[[i]]]
            upDEGs[[names(bottomDEG_idx)[i]]] = rev(rownames(mat)[bottomDEG_idx[[i]]])
        }
        return(list(downDEGs = downDEGs, 
                    upDEGs = upDEGs))
    }

}


#' Run Hierarchical clustering for a matrix
#' 
#' A wrapper for different hierarchical clustering methods to be applied to the within-cell-type
#' cross-conditions Z-score matrix (input). Highly similar conditions (columns) will be grouped together 
#' given the DE Z-scores of rows (genes). 
#' 
#' @export
#' @importFrom protoclust protoclust
#' @importFrom protoclust protocut
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

DEhclust = function(mat = NULL, 
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
        # protoclust.installed <- PackageCheck("protoclust", error = FALSE)
        # if (!protoclust.installed[1]) {
        #     stop("Please install the protoclust package to use MinMax hierarchical clustering method.", 
        #          "\nThis can be accomplished with the following command: ", 
        #          "\n----------------------------------------", "\ninstall.packages('protoclust')", 
        #          "\n----------------------------------------", call. = FALSE)
        # }
        # 
        mat.tree = protoclust::protoclust(mat.dist)
        mat.cut <- protoclust::protocut(mat.tree, h = dist_thres)
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
    return(list(cluster_assignment = cluster_list, 
                full_obj = mat.cut,
                mat = mat))
    
}


#' Run PCApermtest and get significant genes from DEhclust
#' 
#' This function will use the output from the DEhclust() and get the necessary elements for PCApermtest():
#' For each cluster of columns being identified, this function will create a truncated sub-matrix given 
#' the original Z-score matrix. The sub-matrix will only contains the selected columns, and they will be input
#' to the PCApermtest() and get_sig_genes() to get the gene signatures for this cluster. This process will
#' be repeated for each cluster. 
#' 
#' @export
#' @param obj The results object produced by DEhclust() function. 
#' @inheritParams PCApermtest
#' @inheritParams get_sig_genes
#' 
#' @return return a list of vectors, and each vector contains the signature genes identified for each cluster.
#' 

get_sig_genes_DEhclust = function(obj = NULL, 
                                  k = 1, 
                                  var_prop = NULL, 
                                  center = T, scale = T, 
                                  num_iter = 200, 
                                  row_filtering_pval = 0.05, 
                                  var_prop_total = NULL,
                                  perm_pval_thres = 0.05, 
                                  ori_pval_thres = 1.666667e-06, 
                                  cor_threshold = 0.2, 
                                  collapse = T,
                                  ...){
    # check 
    if(length(obj$cluster_assignment) == 0){
        stop("The input obj is not valid. Please check.")
    }
    
    # first run PCApermtest using the output from DEhclust() for each of the cluster 
    res_list = list()
    for(i in 1:length(obj$cluster_assignment)){
        # get the members (column names) for cluster i
        cluster_idx = obj$cluster_assignment[[i]]
        
        # run PCApermtest for them
        perm_res = PCApermtest(mat = obj$mat[, cluster_idx], 
                               k = k, 
                               var_prop = var_prop,
                               var_prop_total = var_prop_total,
                               center = center, scale = scale, 
                               num_iter = num_iter, 
                               row_filtering_pval = row_filtering_pval)
        # get the sig genes given the PCApermtest object
        sig_genes = get_sig_genes(perm_obj = perm_res, 
                                  k = k, 
                                  var_prop = var_prop, 
                                  var_prop_total = var_prop_total,
                                  perm_pval_thres = perm_pval_thres, 
                                  ori_pval_thres = ori_pval_thres, 
                                  cor_threshold = cor_threshold, 
                                  collapse = collapse)
        res_list[[paste0(cluster_idx, collapse = "_")]] = list(sig_genes = sig_genes, 
                                                               perm_res = perm_res)
    }
    return(res_list)
    
}


#' Run MultiCCA for a list of matrices
#' 
#' A function to perform MultiCCA analysis (main function imported from package "PMA", 
#' see PMID 19377034 for details of the algorithm) that takes in a list of multiple 
#' Z-score matrices to find the canonical variates (CVs) that maximize the cross-matrices 
#' correlation. The MultiCCA process is modified so that it is not completed in one 
#' run when multiple rounds of CVs are desired. Instead, after each MultiCCA run, we will
#' identify the columns (samples) from each matrix that highly correlate with the CV of that 
#' matrix, and extract + remove them from the matrix. The next MultiCCA is performed to the 
#' list of such "filtered" matrices. This process is repeated until the desired number of 
#' runs is reached (set by users) or the CVs across the matrices have very low correlation 
#' coefficients. 
#' 
#' 
#' @importFrom PMA MultiCCA
#' @export
#' @param mat_list the list of >= 2 DE Z-score matrices for the multiCCA analysis. Each matrix
#' should have the same named rows, but can have different number of columns (samples). 
#' @param penalty to indicate if penalty should be applied during the multiCCA process. When set
#' to "FALSE", no penalty will be applied so that the results are not forced to be "sparse". However,
#' if a single value or a vector of k values (k should be the same as the number of matrices in the
#' input list), then L1 penalty will be applied to each matrix to force the output CVs to be sparse. 
#' See MultiCCA() from the PMA package for details. 
#' @param standardize a boolen value to indicate whether to standardize each column before running 
#' the MultiCCA. Default is FALSE.
#' @param max_k the maximum number of MultiCCA runs. MultiCCA will be repeated until this number is 
#' reached or the CVs across the matrices have very low correlation coefficients.
#' 
#' @param cor_num for each column of a matrix, the number of CVs that it needs to be significantly correlated 
#' with to be selected as a member of the program. Default is "all", meaning it needs to be significantly 
#' correlated with all other CVs. Alternatively, users may input a integer (e.g., 2, 3, ...).
#' 
#' @param mean_cor_thres During the multiCCA analysis, if any of the input matrix has significantly 
#' low(er) correlation with other matrices, it will impact the MultiCCA process. This is the threshold
#' to remove such matrices. If (the CV of) any matrix has a mean correlation coefficient <= 0.2, it will 
#' be removed and the MultiCCA will be repeated for this round (k-th). Note that this matrix will be appended 
#' back to the list for the next round (k+1-th) of MultiCCA, and the same filtering will be repeated. 
#' Default is 0.2. Set it to 0 to avoid such filtering. Vector is also accepted and will be sequentially used 
#' for each iteration of MultiCCA.
#' 
#' 
#' @seealso [MultiCCA()]
#' @return a list of MultiCCA results for each program identified.

DEmultiCCA = function(mat_list = NULL, 
                      penalty = FALSE, 
                      standardize = FALSE, 
                      max_k = 5, 
                      cor_number = "all",
                      mean_cor_thres = 0.2, 
                      flag_cor_num = T,
                      flag_loose = F,
                      pval_thres = 0.05,
                      cor_coef_thres = 0.6,
                      cor_coef_mean_thres = 0.3,
                      ...) {
    # first, check if PMA is installed. 
    PMA.installed <- PackageCheck("PMA", error = FALSE)
    if (!PMA.installed[1]) {
        stop("Please install the PMA package to use the actual MultiCCA method.", 
             "\nThis can be accomplished with the following command: ", 
             "\n----------------------------------------", "\ninstall.packages('PMA')", 
             "\n----------------------------------------", call. = FALSE)
    }
    
    # check if each of the matrix in the mat_list is named
    if(is.null(names(mat_list))){
        stop("The elements in the mat_list must be named. Please check.")
    }
    
    # check if each matrix in the mat_list has rownames and colnames
    # and check if the rownames of each matrix are the same 
    flag_same_rownames <- TRUE
    for(i in 1:length(mat_list)){
        if(is.null(rownames(mat_list[[i]])) | is.null(colnames(mat_list[[i]]))){
            stop(paste("The", i, "-th matrix in mat_list does not have proper rownames/colnames. Please check."))
        }
        if(i == 1){
            first_matrix_row_names <- rownames(mat_list[[i]])
        } else {
            if (! identical(rownames(mat_list[[i]]), first_matrix_row_names)) {
                same_row_names <- FALSE
                stop(paste("The", i, "-th matrix seems to have different rownames from others. Please check!"))
            }
        }
    }
    
    ##########################################
    #   empty lists to keep track of the multiCCA process
    rm_X = list()      # keep record of the matrix that has mean_cor <= mean_cor_thres
    single_X = list()  # keep record of the matrix that has only 1 column
    res = list()  # the list to store the output from each round of MultiCCA
    
    # 
    cor_coef_thres_list = cor_coef_thres
    cor_coef_mean_thres_list = cor_coef_mean_thres
    
    # run through max_k rounds of MultiCCA
    for(k in 1:max_k){
        # 
        if(length(cor_coef_thres_list) <= k){
            cor_coef_thres = tail(cor_coef_thres_list, 1)
        } else {
            cor_coef_thres = cor_coef_thres_list[k]
        }
        
        # 
        if(length(cor_coef_mean_thres_list) <= k){
            cor_coef_mean_thres = tail(cor_coef_mean_thres_list, 1)
        } else {
            cor_coef_mean_thres = cor_coef_mean_thres_list[k]
        }
        
        ###################
        # this is a flag 
        flag_rm_any_celltype = T
        
        # check if there is any removed matrix from last round 
        # that needs to be appended back to the list
        if(k == 1){
            X2 = mat_list
        } else if(length(rm_X) != 0){
            X2 = c(X2, rm_X)
            rm_X = list()
        } 
        
        # check if the input list of matrices of each round is legit for MultiCCA. 
        # any matrix has mean_cor (across cor_coefs with other matrices) <= mean_cor_thres 
        # will be removed and the MultiCCA will be re-done for this round.
        while (flag_rm_any_celltype == T) {
            # check if there is any matrix only has 1 column. Remove them 
            for(CELLTYPE in names(X2)){
                if(ncol(X2[[CELLTYPE]]) == 1){
                    single_X[[CELLTYPE]] = X2[[CELLTYPE]]
                    X2[[CELLTYPE]] = NULL
                } else if (ncol(X2[[CELLTYPE]]) < 1){
                    X2[[CELLTYPE]] = NULL
                }
            }        
            
            # the names of each matrix
            celltype_list = names(X2)
            
            # set the penalty vector for the MultiCCA (penalty controls the sparsity of the CVs)
            if (penalty == F){
                penalties = sapply(X2, FUN = function(x) sqrt(ncol(x))*1 )
            } else if (length(penalty) == 1) {
                penalties = rep(penalty, length(X2))
            } else {
                stop("Currently we do not support MultiCCA with multiple different penalty values. sorry.")
            }
            
            # check if X2 is empty
            if(length(X2) == 0){
                print(paste("Stopping since the mat_list is now empty."))
                return(list(program_assignment = res, mat_list = mat_list))
            }
            
            # run MultiCCA() from PMA package. 
            out <- PMA::MultiCCA(X2, type=rep("standard",length(X2)),
                            penalty= penalties,niter = 30,
                            ncomponents=1, trace = F, standardize = standardize)
            
            # add col-/row-names to the CVs
            names(out$ws) <- names(X2)
            for(i in names(X2)){
                colnames(out$ws[[i]]) <- paste0("MCP", 1)
                rownames(out$ws[[i]]) <- colnames(X2[[i]])
            }
            
            # get the columns with non-zero loadings if penalty is not set to 1 (sparsity is forced)
            MCP_PRTB_list = list()
            for(CELLTYPE in celltype_list){
                idx = which(out$ws[[CELLTYPE]][,1] != 0)
                MCP_PRTB_list[[CELLTYPE]] = rownames(out$ws[[CELLTYPE]])[idx]
            }
            
            # calculate CCA canonical variates for each matrix
            variates_list = list()
            for(CELLTYPE in celltype_list){
                if(standardize == F){
                    PRTB_variates = as.matrix(X2[[CELLTYPE]]) %*% out$ws[[CELLTYPE]][,1]
                } else {
                    PRTB_variates = scale(as.matrix(X2[[CELLTYPE]])) %*% out$ws[[CELLTYPE]][,1]
                }
                variates_list[[CELLTYPE]] = PRTB_variates
            }
            
            # if single_X is not empty, paste it back to variates
            if(length(single_X) != 0){
                for(CELLTYPE in names(single_X)){
                    variates_list[[CELLTYPE]] = single_X[[CELLTYPE]]
                    X2[[CELLTYPE]] = single_X[[CELLTYPE]]
                    MCP_PRTB_list[[CELLTYPE]] = colnames(single_X[[CELLTYPE]])
                    single_X[[CELLTYPE]] = NULL
                }
            }
            
            # calculate cor_mat and its row_mean (without diagonal elements)
            cor_mat = cor(Reduce(cbind, variates_list), use = "complete.obs")
            # print(cor_mat)
            cor_mat = abs(cor_mat)
            mean_cor = vector()
            for(i in 1:ncol(cor_mat)){
                mean_cor = c(mean_cor, mean(cor_mat[i, -i]))
            }
            
            # now check if any of the mean_cor is <= mean_cor_thres; if so, 
            # remove that matrix from the mat_list and re-do the MultiCCA. 
            if(any(mean_cor < mean_cor_thres)){
                idx_mat_rm = which.min(mean_cor)
                rm_X[[names(X2)[idx_mat_rm]]] = X2[[idx_mat_rm]]
                X2[[idx_mat_rm]] = NULL
            } else {
                flag_rm_any_celltype = F
            }
            
            # stop if X2 only has 1 matrix left 
            if(length(X2) <= 1){
                print(paste("During the", k, "-th MultiCCA, the process has no matrix with mean_cor >=", 
                            mean_cor_thres, "\nTherefore stopping the process."))
                break()
            }
        }
        
        # this is a check after the break() from above loop: 
        # we will stop the whole MultiCCA process here and 
        # return whatever has been obtained from previous rounds of MultiCCA.
        if(length(X2) <= 1){
            print("Stopping the MultiCCA process due to not enough number ( < 2) of correlated matrices! \n")
            max_k = k-1
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
        # get the names of all the matrices (re-doing this step because we might 
        # have removed some matrices)
        celltype_list = names(X2)
        
        # this is an important step: need to define the number of correlation.
        if (cor_number == "all"){
            cor_num = length(celltype_list)
        } else if (cor_number > length(celltype_list)){
            cor_num = length(celltype_list)
        }
        
        # save the MCP list
        # saveRDS(MCP_PRTB_list, file = paste0("MCP", k, "_sparse_", sparsity, "_Parse_", PATHWAY, "_MCP_PRTB_list.rds"))
        
        # save the CCA variates
        # saveRDS(variates_list, file = paste0("Variates_MCP", k, "_sparse_", sparsity, "_Parse_", PATHWAY, ".rds"))
        
        
        ####################################################################
        #  in the following steps, for each matrix we will perform correlation tests between 
        #  its CV and all its columns. We will extract the columns that is significantly 
        #  correlated with at least cor_num CVs. 
        ####################################################################
        
        # some lists to store the intermediate results  
        mean_cor_coef_list = list()
        cor_coef_mat_list = list()
        pvec_mat_BH_list = list()
        pvec_mat_list = list()
        
        for(CELLTYPE in celltype_list){
            pvec_mat = matrix(NA, nrow = length(variates_list), ncol = ncol(X2[[CELLTYPE]]))
            cor_coef_mat = matrix(NA, nrow = length(variates_list), ncol = ncol(X2[[CELLTYPE]]))
            
            # calculate the correlation between each CV and each column
            for(idx in 1:length(variates_list)){
                ## original part 
                pvec = apply(X = X2[[CELLTYPE]], MARGIN = 2,
                             FUN = function(x, y){
                                 return(cor.test(x, y)$p.value)
                             }, 
                             variates_list[[idx]])
                cor_coef = apply(X = X2[[CELLTYPE]], MARGIN = 2,
                                 FUN = cor, 
                                 variates_list[[idx]])
                pvec_mat[idx, ] = pvec
                cor_coef_mat[idx, ] = cor_coef
                colnames(pvec_mat) = colnames(cor_coef_mat) = colnames(X2[[CELLTYPE]])
                
            }
            # print(cor_coef_mat)
            
            # We will adjust the correlation P-values for multiple testing using BH
            pvec_mat_BH = p.adjust(pvec_mat, method = "BH")
            pvec_mat_BH = matrix(pvec_mat_BH, nrow = length(variates_list))
            
            colnames(pvec_mat_BH) = colnames(pvec_mat)
            
            # Aug 23 added
            cor_coef_mat_list[[CELLTYPE]] = cor_coef_mat
            pvec_mat_BH_list[[CELLTYPE]] = pvec_mat_BH
            pvec_mat_list[[CELLTYPE]] = pvec_mat
            
            # back to original procedure
            mean_cor_coef_list[[CELLTYPE]] = colMeans(cor_coef_mat[, , drop = F])
            
        }
        
        # modified on 2023 Jan 23
        rm_row_idx = vector()
        pvec_mat_list_qc = list()
        
        # print(pvec_mat_list)
        # print(cor_coef_mat_list)
        
        for(idx_mat in 1:length(pvec_mat_list)){
            CELLTYPE = names(pvec_mat_list)[idx_mat]
            tmp = pvec_mat_list[[idx_mat]]
            tmp2 = cor_coef_mat_list[[idx_mat]][idx_mat, ]
            tmp3 = colMeans( abs(cor_coef_mat_list[[idx_mat]][-idx_mat, , drop = F]) )
            # round it (2 digits)
            tmp2 = round(x = tmp2, digits = 2)
            tmp3 = round(x = tmp3, digits = 2)
            # print(tmp3)
            
            # 1. p-value thresholding 
            idx_pval = abs(tmp) <= pval_thres
            # 2. cor_coef thresholding
            idx_cor_coef = abs(tmp2) >= cor_coef_thres
            idx_cor_coef = matrix(data = idx_cor_coef, ncol = length(idx_cor_coef), nrow = nrow(tmp), byrow = T)
            # 3. cor_coef_mean_thres 
            idx_cor_coef_mean = abs(tmp3) >= cor_coef_mean_thres
            idx_cor_coef_mean = matrix(data = idx_cor_coef_mean, ncol = length(idx_cor_coef_mean), nrow = nrow(tmp), byrow = T)
            # generate a neigboring graph
            tmp[idx_pval & idx_cor_coef & idx_cor_coef_mean] = 1
            tmp[!(idx_pval & idx_cor_coef & idx_cor_coef_mean)] = 0
            
            # 3. cor_coef thresholding
            # idx_cor_coef2 = abs(tmp2) >= 0.8
            # idx_cor_coef2 = matrix(data = idx_cor_coef2, ncol = length(idx_cor_coef2), nrow = nrow(tmp), byrow = T)
            # tmp[idx_cor_coef2] = 1
            
            # if 
            if(all(tmp[-idx_mat, ] == 0)){
                rm_row_idx = c(rm_row_idx, idx_mat)
            }
            
            colnames(tmp) = paste0(CELLTYPE, "__", colnames(tmp))
            pvec_mat_list_qc[[CELLTYPE]] = tmp
        }
        # print(pvec_mat_list_qc)
        
        if(length(rm_row_idx) != 0){
            # remove the corresponding matrix 
            for(idx_rm in rm_row_idx){
                pvec_mat_list_qc[[idx_rm]] = NULL
            }
            # remove the rows of those corresponding matrix
            for(idx in 1:length(pvec_mat_list_qc)){
                tmp = pvec_mat_list_qc[[idx]]
                tmp = tmp[-rm_row_idx, , drop = F]
                pvec_mat_list_qc[[idx]] = tmp
            }
            # adjust the cor_num
            if(cor_num > length(pvec_mat_list_qc)){
                cor_num = length(pvec_mat_list_qc)
            }
        }
        
        # merge them into one matrix
        test = Reduce(cbind, pvec_mat_list_qc)

        if(length(test) == 0 | ncol(test) <= 1 | nrow(test) <= 1){
            print("No more correlated columns are detected by MultiCCA (flag = 0). Terminating...")
            max_k = k-1
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
        # remove non-relevant column with no correlated CV
        order_1 = colSums(test)
        order_2 = duplicated(t(test)) | duplicated(t(test), fromLast = T)
        test = test[, order(order_1, order_2, decreasing = T)]
        
        test = test[, colSums(test) > 0, drop = F]
        if(length(test) == 0 | ncol(test) <= 1){
            print("No more correlated columns are detected by MultiCCA (flag = 1). Terminating...")
            max_k = k-1
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
        # modified on 2023 Jan 23: 
        
        # here we will assert if cor_num is satisfied; 
        slct_cor_PRTB = vector()
        if (flag_loose == T) {
            if (cor_num > 2){
                slct_cor_PRTB =  colnames(test[, colSums(test) >= (cor_num - 1)])
            } else {
                slct_cor_PRTB =  colnames(test[, colSums(test) >= cor_num ])
            }
        } else if(flag_loose == F & sum(colSums(test) >= cor_num) > 1){
            slct_cor_PRTB = colnames(test[, colSums(test) >= cor_num])
        } else if (flag_loose == F & flag_cor_num == T){
            # we will gradually reduce cor_num by 1.
            cor_num = cor_num - 1
            while(cor_num > 1 & length(slct_cor_PRTB) <= 1){
                slct_cor_PRTB =  colnames(test[, colSums(test) >= cor_num])
                cor_num = cor_num - 1
            }
        } else {
            print("No more correlated columns are detected by MultiCCA (flag = 2). Terminating...")
            max_k = k-1
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
        if(length(slct_cor_PRTB) <= 1){
            print("No more correlated columns are detected by MultiCCA (flag = 3). Terminating...")
            max_k = k-1
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
        # re-store the celltype and PRTB_names
        slct_cor_PRTB = unlist(strsplit(slct_cor_PRTB, split = "__"))
        slct_cor_PRTB = as.data.frame(matrix(slct_cor_PRTB, byrow = T, ncol=2))
        # 
        cor_PRTB_list = list()
        shared_PRTB_list = list()
        for(CELLTYPE in unique(slct_cor_PRTB$V1)){
            cor_PRTB_list[[CELLTYPE]] = slct_cor_PRTB[slct_cor_PRTB$V1 == CELLTYPE, "V2"]
            shared_PRTB_list[[CELLTYPE]] = intersect(cor_PRTB_list[[CELLTYPE]], MCP_PRTB_list[[CELLTYPE]])
        }
        
        # 
        # saveRDS(mean_cor_coef_list, file = paste0("MCP", k, "_sparse_", sparsity, "_Parse_", PATHWAY, "_mean_cor_coef_list.rds"))
        # saveRDS(cor_PRTB_list, file = paste0("MCP", k, "_sparse_", sparsity, "_Parse_", PATHWAY, "_cor_PRTB_list.rds"))
        # saveRDS(shared_PRTB_list, file = paste0("MCP", k, "_sparse_", sparsity, "_Parse_", PATHWAY, "_shared_PRTB_list.rds"))
        
        res[[paste0("Program", k)]] = list(mean_cor_coef_list = mean_cor_coef_list, 
                                           cor_PRTB_list = cor_PRTB_list, 
                                           shared_PRTB_list = shared_PRTB_list)
        
        # remove those PRTBs from mat_list 
        for(CELLTYPE in celltype_list){
            if(length(cor_PRTB_list[[CELLTYPE]]) != 0 & length(X2[[CELLTYPE]]) != 0){
                rm_idx = which(colnames(X2[[CELLTYPE]]) %in% cor_PRTB_list[[CELLTYPE]])
                X2[[CELLTYPE]] = X2[[CELLTYPE]][, -rm_idx, drop=F]
            }
            # 
            if(length(X2[[CELLTYPE]]) == 0){
                X2[[CELLTYPE]] = NULL
                print(paste("\nCelltype", CELLTYPE, "is now empty in mat_list! \n"))
            }
        }
        
        # check if X2 is empty
        if(length(X2) == 0){
            print(paste("Stopping since the mat_list is now empty."))
            return(list(program_assignment = res, mat_list = mat_list))
        }
        
    }
    return(list(program_assignment = res, mat_list = mat_list))
}




#' Run PCApermtest and get significant genes from DEmultiCCA
#' 
#' This function will use the output from the DEmultiCCA() and get the neccessary elements for PCApermtest and 
#' get the gene signatures for each perturbation program that DEmultiCCA() identifies. It works in a similar way 
#' as get_sig_genes_DEhclust() . 
#' 
#' 
#' @export
#' 
#' @param obj The results object produced by DEmultiCCA() function. 
#' @inheritParams PCApermtest
#' @inheritParams get_sig_genes
#' 
#' @return return a list of vectors, and each vector contains the signature genes identified for each MultiCCA program.
#' 

get_sig_genes_DEmultiCCA = function(obj = NULL, 
                                    k = 1, 
                                    var_prop = NULL, 
                                    center = T, scale = T, 
                                    num_iter = 200, 
                                    row_filtering_pval = 0.05, 
                                    var_prop_total = NULL,
                                    perm_pval_thres = 0.05, 
                                    ori_pval_thres = 1.666667e-06, 
                                    cor_threshold = 0.2, 
                                    collapse = T,
                                    ...){
    # check if the input is valid:
    if(length(obj$program_assignment) == 0){
        stop("The input obj is not valid. Please check.")
    }
    
    # first run PCApermtest using the output from DEmultiCCA() for each of the cluster 
    res_list = list()
    
    for(i in 1:length(obj$program_assignment)){
        # get the members for program i
        cluster_info = obj$program_assignment[[i]]$shared_PRTB_list
        
        # loop through all the cell type and get the corresponding columns
        cluster_mat = NULL
        for (CELLTYPE in names(cluster_info)){
            tmp_mat = obj$mat_list[[CELLTYPE]]
            tmp_mat = tmp_mat[ ,cluster_info[[CELLTYPE]], drop = F]
            colnames(tmp_mat) = paste0(CELLTYPE, "__", colnames(tmp_mat))
            if(is.null(cluster_mat)){
                cluster_mat = tmp_mat
            } else {
                cluster_mat = cbind(cluster_mat, tmp_mat)
            }
            rm(tmp_mat)
        }
        
        # run PCApermtest for them
        perm_res = PCApermtest(mat = cluster_mat, 
                               k = k, 
                               var_prop = var_prop,
                               var_prop_total = var_prop_total,
                               center = center, scale = scale, 
                               num_iter = num_iter, 
                               row_filtering_pval = row_filtering_pval)
        # get the sig genes given the PCApermtest object
        sig_genes = get_sig_genes(perm_obj = perm_res, 
                                  k = k, 
                                  var_prop = var_prop, 
                                  var_prop_total = var_prop_total,
                                  perm_pval_thres = perm_pval_thres, 
                                  ori_pval_thres = ori_pval_thres, 
                                  cor_threshold = cor_threshold, 
                                  collapse = collapse)
        res_list[[names(obj$program_assignment)[i]]] = list(sig_genes = sig_genes, 
                                                            perm_res = perm_res)
    }
    return(res_list)
}









