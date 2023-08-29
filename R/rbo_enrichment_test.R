#' Functions for a new gene-set enrichment test based on the 
#' RBO (rank biased overlap) calculation with extropolation (Webber et al., 2010).
#' The core functions of rbo() calculation was modified from the "gespeR" package (original author: Fabian Schmich). 
#' We modified it to accomodate our package and data type. We also developed a permutation scheme for 
#' RBO to allow for p-value calculations. 


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


#' 
#'
#' Function to perform enrichment test based on rbo(). 
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
#' @return a data.frame consists of rbo measurement between the inptu gene list and all the GO terms, 
#' as well as the P-values based on permutaion 
#' 

rbo_enrich_test <- function(input_list, 
                            go_term_db, 
                            p, 
                            n_iter = 500, 
                            k=floor(max(length(list1), length(list2))/2), 
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
    return(res_dat)
}


# 
gs = function(d, p){
    (1-p)*p^(d-1)
}

gs_seq = function(d, p){
    return(gs(1:d, p))
}


