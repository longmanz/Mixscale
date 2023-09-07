#'
NULL


#' Internal Function to Fit a Gamma-Poisson GLM
#'
#' @inheritParams glmGamPoi::glm_gp
#' @inheritParams glmGamPoi::overdispersion_mle
#' @import glmGamPoi
#' 
#' @param Y any matrix-like object (e.g. `matrix()`, `DelayedArray()`, `HDF5Matrix()`) with
#'   one column per sample and row per gene.
#'
#' @return a list with four elements
#'  * `Beta` the coefficient matrix
#'  * `overdispersion` the vector with the estimated overdispersions
#'  * `Mu` a matrix with the corresponding means for each gene
#'     and sample
#'  * `size_factors` a vector with the size factor for each
#'    sample
#'  * `ridge_penalty` a vector with the ridge penalty
#'
#' @seealso [glm_gp()] and [overdispersion_mle()]
#' @keywords internal
glm_gp_disp_only <- function(data,
                   design = ~ 1,
                   col_data = NULL,
                   reference_level = NULL,
                   offset = 0,
                   size_factors = c("normed_sum", "deconvolution", "poscounts", "ratio"),
                   overdispersion = TRUE,
                   overdispersion_shrinkage = TRUE,
                   ridge_penalty = 0,
                   do_cox_reid_adjustment = TRUE,
                   subsample = FALSE,
                   on_disk = NULL,
                   use_assay = NULL,
                   verbose = FALSE){
    
    # Validate `data`
    if(inherits(data, "formula")){
        if(length(design) != 2 || design != ~ 1){
            stop("If the first argument is already a formula, the second argument must not be set. Please call this function like this:\n",
                 "'glm_gp(data = mat, design = ~ a + b + c, ...)'", call. = FALSE)
        }
        extr <- glmGamPoi:::extract_data_from_formula(data, col_data, parent.frame())
        data <- extr$data
        design <- extr$design
    }
    if(is.vector(data)){
        data <- matrix(data, nrow = 1)
    }
    data_mat <- glmGamPoi:::handle_data_parameter(data, on_disk = F)
    
    # Convert the formula to a model_matrix
    col_data <- glmGamPoi:::get_col_data(data, col_data)
    des <- glmGamPoi:::handle_design_parameter(design, data, col_data, reference_level)
    
    # Call glm_gp_impl()
    res <- glm_gp_disp_only_impl(data_mat,
                       model_matrix = des$model_matrix,
                       offset = offset,
                       size_factors = size_factors,
                       overdispersion = overdispersion,
                       overdispersion_shrinkage = overdispersion_shrinkage,
                       ridge_penalty = ridge_penalty,
                       do_cox_reid_adjustment = do_cox_reid_adjustment,
                       subsample = subsample,
                       verbose = verbose)
    # Make sure that the output is nice and beautiful
    names(res$overdispersions) <- rownames(data)

    class(res) <- "glmGamPoi"
    res
}



#' Internal Function to Fit a Gamma-Poisson GLM
#'
#' @inheritParams glm_gp
#' @inheritParams overdispersion_mle
#' @param Y any matrix-like object (e.g. `matrix()`, `DelayedArray()`, `HDF5Matrix()`) with
#'   one column per sample and row per gene.
#'
#' @return a list with four elements
#'  * `Beta` the coefficient matrix
#'  * `overdispersion` the vector with the estimated overdispersions
#'  * `Mu` a matrix with the corresponding means for each gene
#'     and sample
#'  * `size_factors` a vector with the size factor for each
#'    sample
#'  * `ridge_penalty` a vector with the ridge penalty
#'
#' @seealso [glm_gp()] and [overdispersion_mle()]
#' @keywords internal
glm_gp_disp_only_impl <- function(Y, model_matrix,
                        offset = 0,
                        size_factors = c("normed_sum", "deconvolution", "poscounts", "ratio"),
                        overdispersion = TRUE,
                        overdispersion_shrinkage = TRUE,
                        ridge_penalty = 0,
                        do_cox_reid_adjustment = TRUE,
                        subsample = FALSE,
                        verbose = FALSE){
    if(is.vector(Y)){
        Y <- matrix(Y, nrow = 1)
    }
    # Error conditions
    stopifnot(is.matrix(Y) || is(Y, "DelayedArray"))
    stopifnot(is.matrix(model_matrix) && nrow(model_matrix) == ncol(Y))
    glmGamPoi:::validate_Y_matrix(Y)
    subsample <- glmGamPoi:::handle_subsample_parameter(Y, subsample)
    ridge_penalty <- glmGamPoi:::handle_ridge_penalty_parameter(ridge_penalty, model_matrix, verbose = verbose)
    
    # Combine offset and size factor
    off_and_sf <- glmGamPoi:::combine_size_factors_and_offset(offset, size_factors, Y, verbose = verbose)
    offset_matrix <- off_and_sf$offset_matrix
    size_factors <- off_and_sf$size_factors
    
    # Check if there distinct groups in model matrix
    # returns NULL if there would be more groups than columns
    # only_intercept_model <- ncol(model_matrix) == 1 && all(model_matrix == 1)
    groups <- glmGamPoi:::get_groups_for_model_matrix(model_matrix)
    if(! is.null(groups) && any(ridge_penalty > 1e-10)){
        # Cannot apply ridge penalty in group-wise optimization
        groups <- NULL
    }
    
    # If no overdispersion, make rough first estimate
    if(isTRUE(overdispersion)){
        if(verbose){ message("Make initial dispersion estimate") }
        disp_init <- glmGamPoi:::estimate_dispersions_roughly(Y, model_matrix, offset_matrix = offset_matrix)
    }else if(isFALSE(overdispersion)){
        disp_init <- rep(0, times = nrow(Y))
    }else if(is.character(overdispersion) && overdispersion == "global"){
        if(verbose){ message("Make initial dispersion estimate") }
        disp_init <- glmGamPoi:::estimate_dispersions_roughly(Y, model_matrix, offset_matrix = offset_matrix)
        disp_init <- rep(median(disp_init), nrow(Y))
    }else{
        stopifnot(is.numeric(overdispersion) && (length(overdispersion) == 1 || length(overdispersion) == nrow(Y)))
        if(length(overdispersion) == 1){
            disp_init <- rep(overdispersion, times = nrow(Y))
        }else{
            disp_init <- overdispersion
        }
    }
    
    
    # Estimate the betas
    if(! is.null(groups)){
        if(verbose){ message("Make initial beta estimate") }
        beta_group_init <- glmGamPoi:::estimate_betas_roughly_group_wise(Y, offset_matrix, groups)
        if(verbose){ message("Estimate beta") }
        beta_res <- glmGamPoi:::estimate_betas_group_wise(Y, offset_matrix = offset_matrix,
                                              dispersions = disp_init, beta_group_init = beta_group_init,
                                              groups = groups, model_matrix = model_matrix)
    }else{
        # Init beta with reasonable values
        if(verbose){ message("Make initial beta estimate") }
        beta_init <- glmGamPoi:::estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix, ridge_penalty = ridge_penalty)
        if(verbose){ message("Estimate beta") }
        beta_res <- glmGamPoi:::estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                                  dispersions = disp_init, beta_mat_init = beta_init, ridge_penalty = ridge_penalty)
    }
    Beta <- beta_res$Beta
    
    # Calculate corresponding predictions
    # Mu <- exp(Beta %*% t(model_matrix) + offset_matrix)
    Mu <- glmGamPoi:::calculate_mu(Beta, model_matrix, offset_matrix)
    
    # Make estimate of over-disperion
    if(isTRUE(overdispersion) || (is.character(overdispersion) && overdispersion == "global")){
        if(verbose){ message("Estimate dispersion") }
        if(isTRUE(overdispersion)){
            disp_est <- overdispersion_mle(Y, Mu, model_matrix = model_matrix,
                                           do_cox_reid_adjustment = do_cox_reid_adjustment,
                                           subsample = subsample, verbose = verbose)$estimate
        }else if(is.character(overdispersion) && overdispersion == "global"){
            disp_est <- overdispersion_mle(Y, Mu, model_matrix = model_matrix,
                                           do_cox_reid_adjustment = do_cox_reid_adjustment,
                                           global_estimate = TRUE,
                                           subsample = subsample, verbose = verbose)$estimate
            disp_est <- rep(disp_est, times = nrow(Y))
        }
        
        if(isTRUE(overdispersion_shrinkage)){
            dispersion_shrinkage <- overdispersion_shrinkage(disp_est, gene_means = DelayedMatrixStats::rowMeans2(Mu),
                                                             df = subsample - ncol(model_matrix),
                                                             ql_disp_trend  = length(disp_est) >= 100,
                                                             npoints = max(0.1 * length(disp_est), 100),
                                                             verbose = verbose)
            disp_latest <- dispersion_shrinkage$dispersion_trend
        }else{
            dispersion_shrinkage <- NULL
            disp_latest <- disp_est
        }
        
        # Estimate the betas again (only necessary if disp_est has changed)
        # if(verbose){ message("Estimate beta again") }
        # if(! is.null(groups)){
        #     beta_res <- estimate_betas_group_wise(Y, offset_matrix = offset_matrix,
        #                                           dispersions = disp_latest, beta_mat_init = Beta,
        #                                           groups = groups, model_matrix = model_matrix)
        # }else{
        #     beta_res <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
        #                                               dispersions = disp_latest, beta_mat_init = Beta, ridge_penalty = ridge_penalty)
        # }
        # Beta <- beta_res$Beta
        # 
        # # Calculate corresponding predictions
        # Mu <- calculate_mu(Beta, model_matrix, offset_matrix)
    }else if(isTRUE(overdispersion_shrinkage) || is.numeric(overdispersion_shrinkage)){
        # Given predefined disp_est shrink them
        disp_est <- disp_init
        dispersion_shrinkage <- overdispersion_shrinkage(disp_est, gene_means = DelayedMatrixStats::rowMeans2(Mu),
                                                         df = subsample - ncol(model_matrix),
                                                         disp_trend = overdispersion_shrinkage, verbose = verbose)
        disp_latest <- dispersion_shrinkage$dispersion_trend
        # if(verbose){ message("Estimate beta again") }
        # if(! is.null(groups)){
        #     beta_res <- estimate_betas_group_wise(Y, offset_matrix = offset_matrix,
        #                                           dispersions = disp_latest, beta_mat_init = Beta,
        #                                           groups = groups, model_matrix = model_matrix)
        # }else{
        #     beta_res <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
        #                                               dispersions = disp_latest, beta_mat_init = Beta, ridge_penalty = ridge_penalty)
        # }
        # Beta <- beta_res$Beta
        # # Calculate corresponding predictions
        # Mu <- calculate_mu(Beta, model_matrix, offset_matrix)
    }else{
        # Use disp_init, because it is already in vector shape
        disp_est <- disp_init
        dispersion_shrinkage <- NULL
    }
    
    
    # Return everything
    list(Beta = Beta,
         overdispersions = disp_est,
         overdispersion_shrinkage_list = dispersion_shrinkage)
}



