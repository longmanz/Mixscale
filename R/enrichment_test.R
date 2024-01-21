#'
NULL

#' Wrapper function for DE and enrichment test
#' 
#' This function provides a wrapper of Seurat::FindMarkers() and Mixscale::fisher_enrich_test(). 
#' Users can input a Seurat object they want to investigate and a list of gene sets they want to 
#' test against, and the wrapper will perform DE tests + Fisher's enrichment test across all the 
#' available cell types. It will then return a list of data frames, containing gene set enrichment 
#' results for each cell type. 
#' 
#' @export
#' @param object a seurat object to perform the DE test and the enrichment test
#' @param plist the pathway gene lists to test the DE genes against
#' @param ct.class Regroup cells into a different identity class prior to performing differential expression. 
#' Default is NULL (so all cells be used simultaneously).
#' @param slct.ct Subset a particular identity class prior to regrouping. Only relevant if group.by is set.
#' @param ident.1 Identity class to define markers for; pass an object of class phylo or 'clustertree' to find markers for a node in a cluster tree; passing 'clustertree' requires BuildClusterTree to have been run
#' @param ident.2 A second identity class for comparison; if NULL, use all other cells for comparison; if an object of class phylo or 'clustertree' is passed to ident.1, must pass a node to find markers for
#' @return a list of data frames containing the gene set enrichment results for each group in "group.by"
#' 

DEenrich <- function(object, 
                     plist = NULL, 
                     ident = NULL, 
                     ident.1 = NULL,
                     ident.2 = NULL, 
                     ct.class = NULL, 
                     slct.ct = NULL,
                     direction = c("up", "down", "both"), 
                     logfc.threshold = 0.25,
                     p.val.cutoff = 0.05, 
                     min.pct = 0.1,
                     assay = NULL, 
                     ...){
    
    slct_celltype = sort(unique(object[[ct.class]][, 1]))
    if(!is.null(ct.class) & is.null(slct_celltype)){
        stop("Please check if your ct.class is correctly specified.")
    }
    
    if(!is.null(slct.ct)){
        slct_celltype = intersect(slct_celltype, slct.ct)
    }
    if(length(slct_celltype) == 0){
        slct_celltype = "con1"
    }
    
    enrich_list = list()
    for(CELLTYPE in slct_celltype){
        if(!is.null(ident)){
            Idents(object) = ident
        }
        # 
        if(is.null(ct.class)){
            object[["new_ident"]] = paste0("con1", "_", Idents(object))
        } else {
            object[["new_ident"]] = paste0(object[[ct.class]][,1], "_", Idents(object))
        }
        ident.1.tmp = paste0(CELLTYPE, "_", ident.1)
        ident.2.tmp = paste0(CELLTYPE, "_", ident.2)
        
        # run DE 
        Idents(object) = "new_ident"
        DE_res = FindMarkers(object, 
                             ident.1 = ident.1.tmp, 
                             ident.2 = ident.2.tmp, 
                             min.pct = min.pct, 
                             logfc.threshold = 0, 
                             group.by = group.by,
                             subset.ident = subset.ident,
                             ...)
        
        # get the top DEGs separately for up and down regulated genes 
        upDEG = rownames(DE_res[DE_res$p_val_adj <= p.val.cutoff & DE_res$avg_log2FC > logfc.threshold, ])
        downDEG = rownames(DE_res[DE_res$p_val_adj <= p.val.cutoff & DE_res$avg_log2FC < -logfc.threshold, ])
        background = rownames(DE_res) # the background gene list
        
        # run enrichment test for the DEGs
        if(length(downDEG) < 5 | direction == "up"){
            enrich_res_down = NULL
        } else {
            enrich_res_down = fisher_enrich_test(input_list = downDEG, 
                                                 background = background, 
                                                 go_term_db = plist)
            enrich_res_down$num_DEG = length(downDEG)
            enrich_res_down$direction_DEG = "downDEG"
            enrich_res_down = enrich_res_down[order(enrich_res_down$Pval), ]
            
        }

        # 
        if(length(upDEG) < 5 | direction == "down"){
            enrich_res_up = NULL
        } else {
            enrich_res_up = fisher_enrich_test(input_list = upDEG, 
                                               background = background, 
                                               go_term_db = plist)
            enrich_res_up$num_DEG = length(upDEG)
            enrich_res_up$direction_DEG = "upDEG"
            enrich_res_up = enrich_res_up[order(enrich_res_up$Pval), ]
        }
        
        # save the results to the list
        enrich_list[[CELLTYPE]] = rbind(enrich_res_up, enrich_res_down)
    }
    
    if(length(enrich_list) == 1){
        return(enrich_list[[1]])
    } else {
        return(enrich_list)
    }
}



#' Rank biased overlap 
#' 
#' A function for a new gene-set enrichment test based on the 
#' RBO (rank biased overlap) calculation with extropolation (Webber et al., 2010).
#' The core functions of rbo() calculation was modified from the "gespeR" package (original author: Fabian Schmich). 
#' We modified it to accomodate our package and data type. We also developed a permutation scheme for 
#' RBO to allow for p-value calculations. 
#' 
#' @author Fabian Schmich ("gespeR" package)
#' @export
#' 
#' @param list1 List 1
#' @param list2 List 2
#' @param p Weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k Evaluation depth for extrapolation
#' @param side Evaluate similarity between the top or the bottom of the ranked lists
#' @param mid Set the mid point to for example only consider positive or negative scores
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return a scaler value measuring the rank biased overlap (rbo)
#' 

rbo <- function(list1, list2, p, k=floor(max(length(list1), length(list2))/2), side=c("top", "bottom"), mid = NULL, uneven.lengths = TRUE) {
    side <- match.arg(side)
    if (!is.numeric(list1) | !is.numeric(list2))
        stop("Input vectors are not numeric.")
    if (is.null(names(list1)) | is.null(names(list2)))
        stop("Input vectors are not named.")
    ids <- switch(side,
                  "top"=list(list1=.select.ids(list1, "top", mid), list2=.select.ids(list2, "top", mid)),
                  "bottom"=list(list1=.select.ids(list1, "bottom", mid), list2=.select.ids(list2, "bottom", mid))
    )
    min(1, .rbo.ext(ids$list1, ids$list2, p, k, uneven.lengths = uneven.lengths))
}


# rbo2 <- function(list1, list2, p, k=floor(max(length(list1), length(list2))/2), side=c("top", "bottom"), mid = NULL, uneven.lengths = TRUE) {
#     side <- match.arg(side)
#     if (!is.numeric(list1) | !is.numeric(list2))
#         stop("Input vectors are not numeric.")
#     if (is.null(names(list1)) | is.null(names(list2)))
#         stop("Input vectors are not named.")
#     ids <- switch(side,
#                   "top"=list(list1=.select.ids(list1, "top", mid), list2=.select.ids(list2, "top", mid)),
#                   "bottom"=list(list1=.select.ids(list1, "bottom", mid), list2=.select.ids(list2, "bottom", mid))
#     )
#     min(1, rbo_ext(ids$list1, ids$list2, p, k, uneven_lengths = uneven.lengths))
# }


#' Select top or bottom names of ranked vector
#' 
#' @author Fabian Schmich ("gespeR" package)
#' @noRd
#' 
#' @param x The ranked list
#' @param side The side to be evaluated ("top" or "bottom" of ranked list)
#' @param mid The mid point to split a list, e.g. to split between positive and negative values choose mid=0
#' @return A vector of selected identifiers
.select.ids <- function(x, side=c("top", "bottom"), mid=NULL) {
    side <- match.arg(side)
    if (side == "top")  {
        x <- sort(x, decreasing=TRUE)
        if (is.null(mid))
            return(names(x))
        else 
            return(names(x)[which(x > mid)])
    } else if (side == "bottom") {
        x <- sort(x, decreasing=FALSE)
        if (is.null(mid)) 
            return(names(x))
        else 
            return(names(x)[which(x < mid)])
    }
}


#' Rank biased overlap formula based on (32) from "A Similarity Measure for Indefinite Rankings" (Webber et al.)
#' 
#' @author Fabian Schmich ("gespeR" package)
#' @noRd
#' 
#' @param x List 1
#' @param y List 2
#' @param p The weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k The evaluation depth
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return The rank biased overlap between x and y
.rbo.ext <- function(x, y, p, k, uneven.lengths = TRUE) {
    if (length(x) <= length(y)) {
        S <- x
        L <- y
    } else {
        S <- y
        L <- x
    }
    l <- min(k, length(L))
    s <- min(k, length(S))
    
    if (uneven.lengths) {
        Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
        ((1-p) / p) *
            ((sum(Xd[seq(1, l)] / seq(1, l) * p^seq(1, l))) +
                 (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
            ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l  
    } else {
        #stopifnot(l == s)
        k <- min(s, k)
        Xd <- sapply(1:k, function(i) length(intersect(x[1:i], y[1:i])))
        Xk <- Xd[k]
        (Xk / k) * p^k + (((1-p)/p) * sum((Xd / seq(1,k)) * p^seq(1,k)))
    }
}


#' Rank biased overlap (RBO) based enrichment test
#'
#' To perform enrichment test based on rank biased overlap and permutation. 
#'
#' @export
#' 
#' @param input_list input gene list from user (a named vector)
#' @param go_term_db a list object of multiple gene-ontology (GO) terms to run enrichment test against
#' @param p Weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param n_iter the number of iteration to perform the permutation to obtain the P-values of the enrichment test 
#' @param k Evaluation depth for extrapolation
#' @param side Evaluate similarity between the top or the bottom of the ranked lists
#' @param mid Set the mid point to for example only consider positive or negative scores
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @param empirical_test a boolen value to tell the function is an empirical test should be performed. If TRUE,
#' the exact empirical proportion of the permutated elements that are greater than the true RBO 
#' is returned as the p-value (high accuracy usually requires a large n_iter, e.g., 1000). If FALSE, then a standard 
#' Z-score test is applied to the RBO based on the mean and standard deviation of all the permuated elements (less accurate 
#' but more efficient. A small n_iter is usually enough (e.g., 100 or 200) to get good approximation compared to 
#' the true empirical test). 
#' 
#' @return a data.frame consists of rbo measurement between the inptu gene list and all the GO terms, 
#' as well as the P-values based on permutation. Please note that the P-values indicate whether the rank of the input gene
#' list and the GO-term gene set are consistent or not. It does NOT indicate if RBO is significantly different from 0. 
#' 

rbo_enrich_test <- function(input_list, 
                            go_term_db, 
                            p, 
                            n_iter = 500, 
                            k=300, 
                            side=c("top", "bottom"), 
                            mid = NULL, 
                            uneven.lengths = TRUE, 
                            empirical_test = FALSE, 
                            seed = 131415926) {
    if(length(input_list) < 5){
        print("The length of the input list is less than 5, stopping analysis...")
        return(NULL)
    }
    
    # 0. each vector in the go_term_db is assumed to be an ordered vector of characters (gene names). 
    #    we need to convert each vector in to a named vector of number (number being the rank of each gene). 
    go_term_db2 = lapply(X = go_term_db, 
                         FUN = function(x) {
                             if(side == "bottom") tmp = 1:length(x)
                             else if (side == "top") tmp = length(x):1
                             names(tmp) = x
                             return(tmp)
                         })
    
    # create a rank vector from input_list
    ori_input_list = input_list
    if(side == "bottom"){
        input_list = 1:length(ori_input_list)
        names(input_list) = ori_input_list
    } else if (side == "top") {
        input_list = length(ori_input_list):1
        names(input_list) = ori_input_list
    }
    
    # 1. calculate the true RBO for the input_list and all the GO terms
    rbo_real = sapply(X = go_term_db2, 
                      FUN = rbo, 
                      list2 = input_list, 
                      p = p, 
                      k = k, 
                      side = side, 
                      mid = mid, 
                      uneven.lengths = uneven.lengths)
    
    # 2. now we need to proceed to the permutation tests 
    set.seed(seed)
    
    # Shuffle the input_list for n_iter times (shuffled matrix has n_iter columns)
    max_d = ifelse(k >= length(input_list), yes = length(input_list), no = k)
    shuffled_matrix <- replicate(n_iter, sample(input_list[1:max_d]))
    rownames(shuffled_matrix) = names(input_list[1:max_d])
    
    # Function to calculate rbo for each go_term against the shuffled matrix
    calculate_rbo <- function(go_term) {
        apply(shuffled_matrix, 
              MARGIN = 2, 
              FUN = rbo, 
              list2 = go_term, 
              p = p, 
              k = k, 
              side = side, 
              mid = mid, 
              uneven.lengths = uneven.lengths)
    }
    
    # Apply the function to each go_term in go_term_db2
    list_perm_vector <- lapply(go_term_db2, calculate_rbo)
    
    # now calculate the P-values based on empirical_test
    if(empirical_test == TRUE){
        calculate_proportion <- function(element, list_vec) {
            mean(list_vec > element)
        }
        p_values <- mapply(calculate_proportion, rbo_real, list_perm_vector)
    } else if (empirical_test == FALSE){
        calculate_proportion <- function(element, list_vec) {
            pnorm(q = element, 
                  mean = mean(list_vec), 
                  sd = sd(list_vec), 
                  lower.tail = F)
        }
        p_values <- mapply(calculate_proportion, rbo_real, list_perm_vector)
    }
    
    # merge the results into one data.frame
    res_dat = data.frame(GO_term = names(rbo_real), 
                         RBO = rbo_real, 
                         Pval = p_values, 
                         n_GO_term = sapply(X = go_term_db2, FUN = length), 
                         n_intersect = sapply(X = go_term_db, FUN = function(x, y) length(intersect(x, y)), y = names(input_list)))
    # 
    res_dat$Pval[res_dat$n_intersect <= 1 | res_dat$RBO <= 0.01 ] = 1
    rownames(res_dat) = NULL
    return(res_dat)
}


#' get the weight of each depth till 'd' given a weight parameter 'p'
#' @noRd
#' @return a numeric vector of the weights from rank depth 1 to d.
gs_seq = function(d, p){
    gs = function(d, p){
        (1-p)*p^(d-1)
    }
    return(gs(1:d, p))
}



#' Standard Fisher's exact test for enrichment analysis
#' 
#' This function will perform the strandard Fisher's exact test between the input gene 
#' list and a series of gene-ontology gene sets (adopted from DAVID GO analysis).
#' 
#' @export
#' 
#' @param input_list the input gene list 
#' @param background the background gene list (usually the expressed genes where the 
#' input gene list is generate from, ).
#' @param go_term_db a list of gene-lists (GO term). It should be a list contain multiple named vector, 
#' and each vector should be a vector of multiple marker/signature genes for some biological pathway/process.
#' @param list_gene A Boolen value to indicate if the overlapping genes between the input gene list and 
#' the GO-term should be output as well.
#' @param EASE A Boolen value to indicate if the EASE correction should be applied (see 
#' https://david.ncifcrf.gov/helps/functional_annotation.html). This is useful to mitigate the 
#' small-sample inflation when the input gene list is short (e.g., < 10).
#' 
#' @return a data frame contains the enrichment test results. Each row contains the P-value and enrichment odds
#' ratio calculated from a Fisher's exact test for one GO-term in the go_term_db. 

fisher_enrich_test = function(input_list = NULL, 
                   background = NULL, 
                   go_term_db = NULL, 
                   list_gene = F, 
                   EASE = F){
    PT = length(background)
    
    # if go_term_db is a list of lists, Reduce it down to a list of vector (remove all the intermediate layers)
    while(class(go_term_db[[1]]) == "list"){
        go_term_db = Reduce(c, go_term_db)
    }
    
    if(list_gene == F){
        dat = matrix(nrow = length(go_term_db), ncol = 6)
        i = 1
        for(GO_TERM in names(go_term_db)){
            PH = length(intersect(go_term_db[[GO_TERM]], background))
            LT = length(input_list)
            
            # LH_list = intersect(input_list, go_term_db[[GO_TERM]])
            LH = length(intersect(input_list, go_term_db[[GO_TERM]]))
            
            # the Fisher exact test with EASE correction
            dat2 = matrix(data = c(LH-1, PH-LH+1, LT-LH, PT-LT-(PH-LH)), byrow = T, ncol = 2)
            # print(dat2)
            if(LH < 1 | EASE == T){
                dat2 = matrix(data = c(LH, PH-LH, LT-LH, PT-LT-(PH-LH)), byrow = T, ncol = 2)
            }
            # 
            res = fisher.test(dat2, alternative = "greater")
            # 
            dat[i, 1] = GO_TERM
            dat[i, 2] = res$estimate
            dat[i, 3] = res$p.value
            dat[i, 4] = -log(res$p.value)*res$estimate
            dat[i, 5] = LH
            dat[i, 6] = PH
            # 
            i = i + 1
        }
        dat = as.data.frame(dat)
        dat = dat[complete.cases(dat), ]
        names(dat) = c("GO_term", "OR", "Pval", "combined_score", "num_LH", "num_PH")
        
    } else {
        dat = matrix(nrow = length(go_term_db), ncol = 7)
        i = 1
        for(GO_TERM in names(go_term_db)){
            PH = length(intersect(go_term_db[[GO_TERM]], background))
            LT = length(input_list)
            
            LH_list = intersect(input_list, go_term_db[[GO_TERM]])
            LH = length(LH_list)
            
            # the Fisher exact test with EASE correction
            dat2 = matrix(data = c(LH-1, PH-LH+1, LT-LH, PT-LT-(PH-LH)), byrow = T, ncol = 2)
            # print(dat2)
            if(LH < 1){
                dat2 = matrix(data = c(LH, PH-LH, LT-LH, PT-LT-(PH-LH)), byrow = T, ncol = 2)
            }
            # 
            res = fisher.test(dat2, alternative = "greater")
            # 
            dat[i, 1] = GO_TERM
            dat[i, 2] = res$estimate
            dat[i, 3] = res$p.value
            dat[i, 4] = -log(res$p.value)*res$estimate
            dat[i, 5] = LH
            dat[i, 6] = PH
            dat[i, 7] = paste0(LH_list, collapse = ";")
            
            # 
            i = i + 1
        }
        dat = as.data.frame(dat)
        dat = dat[complete.cases(dat), ]
        names(dat) = c("GO_term", "OR", "Pval", "combined_score", "num_LH", "num_PH", "overlap_gene")
        
    }
    
    dat$OR = as.numeric(dat$OR)
    dat$Pval = as.numeric(dat$Pval)
    dat$num_LH = as.integer(dat$num_LH)
    dat$num_PH = as.integer(dat$num_PH)
    dat$n_GO_term = sapply(X = go_term_db, FUN = length)
    
    return(dat)
}

















