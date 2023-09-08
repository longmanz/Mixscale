#' @include get_fold_change.R
#' @include glm_gp_disp_only.R
#'
NULL

#' Scoring-based DE test
#'
#' A function to perform differential expression (DE) tests based on the perturbation scores from the 
#' PRTBScoring() function. It is a multivariate negative binomial based model that incorporates both the heterogeneity 
#' of perturbation strength in each cell, as well as their cell type background. 
#'
#' @inheritParams Seurat::FindMarkers
#' @import Seurat
#' @import SeuratObject
#' @import glmGamPoi
#' @importFrom Matrix rowMeans
#' 
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene labels.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param PRTB_list provide a vector of PRTBs that the DE tests are restricted 
#' to. Default is NULL (DE tests will be performed for all available PRTBs).
#' @param new.class.name Name of mixscape classification to be stored in
#' metadata.
#' @param min.de.genes Required number of genes that are differentially
#' expressed for method to separate perturbed and non-perturbed cells.
#' @param min.cells Minimum number of cells in target gene class. If fewer than
#' this many cells are assigned to a target gene class during classification,
#' all are assigned NP.
#' @param split.by metadata column with experimental condition/cell type
#' classification information. This is meant to be used to account for cases a
#' perturbation is condition/cell type -specific.
#' @param logfc.threshold the log-fold-change threshold to select genes for DE
#' test. Genes with log-fold-change larger than this value will be selected for DE test. 
#' Note that if split.by is set and more than 1 split.by group exists, this 
#' logfc.threashold will be applied to each group. Default is 0.25.
#' @param total_ct_labels metadata column for the total RNA counts of each cell
#' 
#' @return list of DE results 
#' @export
#' @concept perturbation_scoring

scoringDE = function (object, assay = "RNA", slot = "data", labels = "gene", 
                      nt.class.name = "NT", verbose = FALSE, 
                      PRTB_list = NULL,
                      split.by = "cell_type",
                      total_ct_labels = "nCount_RNA",  
                      logfc.threshold = 0.2, 
                      pseudocount.use = 0.1, 
                      base = 2,
                      min.pct = 0.1, 
                      min.cells = 10) 
{
    # 
    print("Running scoring DE test!\n")
    
    # check if the required package is installed
    glmGamPoi.installed <- SeuratObject::PackageCheck("glmGamPoi", error = FALSE)
    if (!glmGamPoi.installed[1]) {
        stop("Please install the glmGamPoi package to use scoringDE", 
             "\nThis can be accomplished with the following command: ", 
             "\n----------------------------------------", "\nBiocManager::install('glmGamPoi')", 
             "\n----------------------------------------", call. = FALSE)
    }
    
    # 
    assay <- assay %||% DefaultAssay(object = object)
    if (is.null(x = labels)) {
        stop("Please specify target gene class metadata name")
    }
    
    # get the PRTBs with weights
    prtb_score <- Tool(object = object, slot = "PRTBScoring")
    wt_PRTB_list = sort(names(prtb_score))
    
    # get the full list of PRTBs 
    all_PRTB_list = sort(unique(object[[labels]][,1]))
    all_PRTB_list = all_PRTB_list[all_PRTB_list != nt.class.name]
    wt_PRTB_list = wt_PRTB_list[wt_PRTB_list %in% all_PRTB_list]
    DefaultAssay(object) = "RNA"
    
    if(!is.null(PRTB_list)){
        wt_PRTB_list = intersect(wt_PRTB_list, PRTB_list)
        all_PRTB_list = intersect(all_PRTB_list, PRTB_list)
        # check if there is any PRTB left for DE test
        if(length(all_PRTB_list) == 0){
            stop("There is no perturbation left for DE tests. Please check the PRTB_list you provided.")
        }
    }
    
    # 
    # count_data = GetAssayData(object = object[['RNA']], slot = "counts")   # get the raw count data (required by neg binom)
    # count_data_std = GetAssayData(object = object[['RNA']], slot = "data")   # get the std data (required by log-fold change test)
    
    # get the overall covariate matrix from the meta.data
    mat_B = data.frame(cell_label = colnames(object), 
                       nCount_RNA = object[[total_ct_labels]][,1], 
                       cell_type = as.factor(object[[split.by]][,1]), 
                       gene = object[[labels]][,1] )
    
    # the list to store all the DE results for every PRTB
    all_res = list()
    
    # perform scoring DE test for each PRTB
    for(PRTB in all_PRTB_list){
        # first we need to extract the perturbation score for this PRTB (if it has the scores)
        mat_A = data.frame()
        
        if(PRTB %in% wt_PRTB_list){
            celltype_list = names(prtb_score[[PRTB]])
            if(verbose){
                print(paste0(PRTB, " has good scores. will use actual score as weights."))
            }
            for(celltype in celltype_list){
                tmp = prtb_score[[PRTB]][[celltype]]
                
                # get the idx for NT cells and PRTBed cells
                idx_NT = which(tmp$gene == "NT")
                idx_gene = which(tmp$gene == PRTB)
                
                #  1. calculate the overall weights
                # calculate the mean and sd of the PRTB score for the NT cells
                mean_NT = mean(tmp$pvec[idx_NT], na.rm = T)
                sd_NT = sd(tmp$pvec[idx_NT], na.rm = T)
                # standardize the PRTB scores for the PRTBed cells based on the mean and SD from those of the NT cells
                std_weight_gene = (tmp$pvec[idx_gene] - mean_NT)/sd_NT
                # convert those negative standardised PRTB score to 0
                std_weight_gene[which(std_weight_gene < 0)] = 0
                
                # create a new column called "weight" in the tmp dataframe. 
                tmp$weight = 0 
                tmp$weight[idx_gene] = std_weight_gene
                
                # now calculate the LOO weights
                for(idx_col in 3:(ncol(tmp) - 1)){
                    rm(mean_NT, sd_NT, std_weight_gene)
                    
                    # calculate the mean and sd of the PRTB score for the NT cells
                    mean_NT = mean(tmp[idx_NT, idx_col], na.rm = T)
                    sd_NT = sd(tmp[idx_NT, idx_col], na.rm = T)
                    
                    # standardize the PRTB scores for the PRTBed cells based on the mean and SD from those of the NT cells
                    std_weight_gene = (tmp[idx_gene, idx_col] - mean_NT)/sd_NT
                    # convert those negative standardised PRTB score to 0
                    std_weight_gene[which(std_weight_gene < 0)] = 0
                    
                    # create a new column called "weight" in the tmp dataframe. 
                    tmp[paste0("weight_", colnames(tmp)[idx_col])] = 0 
                    tmp[idx_gene, paste0("weight_", colnames(tmp)[idx_col])] = std_weight_gene
                    
                }
                
                tmp$cell_label = row.names(tmp)
                tmp = tmp[, c("cell_label", "pvec", "gene", grep(pattern = "weight", x = colnames(tmp), value = T))]
                
                mat_A = rbind(mat_A, tmp)
                rm(tmp)
            }
            
            ##############################################
            # . merge the mat_A and mat_B , 
            mat_all = merge(mat_A, mat_B, by = c("cell_label", "gene"))
            rm(mat_A)
            mat_all$log_ct = log10(mat_all$nCount_RNA)
            
            # DE Flag:
            DE_FLAG = "weighted"
        } else {
            celltype_list = names(prtb_score[[1]])
            
            print(paste0(PRTB, " does not have score. Will use 0/1 for coding."))
            # 
            tmp = mat_B[mat_B$gene %in% c(PRTB, "NT"), ]
            tmp$weight = 0
            tmp[tmp$gene == PRTB, "weight"] = 1
            
            mat_all = tmp 
            rm(tmp)
            mat_all$log_ct = log10(mat_all$nCount_RNA)
            mat_all = mat_all[order(mat_all$weight), ]
            rownames(mat_all) = 1:nrow(mat_all)
            
            # DE Flag:
            DE_FLAG = "standard"
        }
        
        ##############################################
        # get idx for the cells of this PRTB
        idx = match(mat_all$cell_label, colnames(object))
        
        # idx_NT = match(mat_all$cell_label[mat_all$gene == "NT"], colnames(count_data2))
        # idx_P = match(mat_all$cell_label[mat_all$gene != "NT"], colnames(count_data2))
        
        count_data2 = GetAssayData(object = object[['RNA']], slot = "counts")[, idx]
        count_data_std2 = GetAssayData(object = object[['RNA']], slot = "data")[, idx]
        
        # do not do fold-change check, only do var and min.pct check
        # idx_for_DE = which(apply(count_data_std2, MARGIN = 1, FUN = get_idx, idx_P = idx_P, idx_NT = idx_NT, logfc.threshold = logfc.threshold, norm.method = 'log.norm'))
        # overall_FC = FoldChange_new(object = GetAssayData(object = object, slot = "data"),
        #                                 cells.1 = colnames(object)[idx_NT],
        #                                 cells.2 = colnames(object)[idx_P],
        #                                 mean.fxn = function(x) log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base),
        #                                 fc.name = "avg_log2FC",
        #                                 features = rownames(x = object) ) 
        # 
        # idx_for_DE = (overall_FC$avg_log2FC != 0) & 
        #     (overall_FC$pct.1 >= min.pct | overall_FC$pct.2 >= min.pct) & 
        #     (overall_FC$min.cell.1 >= min.cells | overall_FC$min.cell.2 >= min.cells)
        
        # here we will calculate the logfc within each cell type. 
        fc_list = list()
        for( idx_i in 1:(length(celltype_list)) ){
            celltype = levels(mat_all$cell_type)[idx_i]
            # get the indices
            idx_NT_celltype = match(mat_all$cell_label[mat_all$cell_type == celltype & 
                                                           mat_all$gene == "NT"], colnames(object))
            idx_P_celltype = match(mat_all$cell_label[mat_all$cell_type == celltype & 
                                                          mat_all$gene != "NT"], colnames(object))
            # get the fold-change and min.pct and min.cell
            fc = FoldChange_new(obj = GetAssayData(object = object, slot = "data"),
                                cells.1 = colnames(object)[idx_NT_celltype],
                                cells.2 = colnames(object)[idx_P_celltype],
                                mean.fxn = function(x) log(x = Matrix::rowMeans(x = expm1(x = x)) + pseudocount.use, base = base),
                                fc.name = "avg_log2FC",
                                features = rownames(x = object) ) 
            # overall filtering (including fold-change)
            fc$status = fc$avg_log2FC >= logfc.threshold & 
                (fc$pct.1 >= min.pct | fc$pct.2 >= min.pct) & 
                (fc$min.cell.1 >= min.cells | fc$min.cell.2 >= min.cells)
            
            # filtering on pct and min.cells (genes will be filtered if these two criteria are not met)
            fc$status2 = (fc$pct.1 >= min.pct | fc$pct.2 >= min.pct) & 
                (fc$min.cell.1 >= min.cells | fc$min.cell.2 >= min.cells)
            
            fc_list[[celltype]] = fc 
            rm(fc)
        }
        
        # loop through all the cell types in the fc_list to get a index for the genes to test 
        for( idx_i in 1:(length(celltype_list)) ){
            celltype = levels(mat_all$cell_type)[idx_i]
            
            if(idx_i == 1){
                idx_list = data.frame(celltype = fc_list[[celltype]]$status)
                idx_list2 = data.frame(celltype = fc_list[[celltype]]$status2)
                names(idx_list) = celltype
                # 
                fc_mat = data.frame(celltype = fc_list[[celltype]]$avg_log2FC)
                names(fc_mat) = paste0("fc_", celltype)
            } else {
                idx_list[[celltype]] = fc_list[[celltype]]$status
                idx_list2[[celltype]] = fc_list[[celltype]]$status2
                
                fc_mat[[paste0("fc_", celltype)]] = fc_list[[celltype]]$avg_log2FC
            }
        }
        # count all the columns and get a single index vector
        idx_for_DE = which(apply(X = idx_list, MARGIN = 1, FUN = any))
        # get the fold-change matrix
        fc_mat[!as.matrix(idx_list2)] = NA
        
        # need to add the PRTB target itself
        idx_PRTB = which(rownames(object) %in% PRTB)
        # merge the two indices
        idx_for_DE = unique(c(idx_PRTB, idx_for_DE))
        
        
        ######################################################
        ###  the actual DE test using glmGamPoi
        ######################################################
        
        if(DE_FLAG == "weighted"){
            print(paste0(PRTB, " is running weighted DE test! "))
            
            ######################################################
            # since we will use Leave-one-out to run the DE, we will need some more filterings
            de_gene_to_rm = gsub(pattern = "weight_", replacement = "", x = grep(pattern = "weight_", x = colnames(mat_all), value = T))
            de_gene_to_rm = unique(c(PRTB, de_gene_to_rm) )
            
            # this is the index for genes using LOO DE (see below for the 2nd DE step)
            idx_LOO_rm = which(row.names(object) %in% de_gene_to_rm)
            # this is the filtered index for genes using standard DE (the 1st DE step)
            idx_for_DE = setdiff(idx_for_DE, idx_LOO_rm)
            
            ##################################
            ##    1.  standard DE test 
            # fit the model
            if(length(celltype_list) > 1){
                fit_rough <- glm_gp_disp_only(data = count_data2[idx_for_DE, ], 
                                    design = ~ 0 + cell_type + log_ct, 
                                    col_data = mat_all, 
                                    size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                    on_disk = FALSE)
                # 
                fit <- glm_gp(data = count_data2[idx_for_DE, ], 
                              design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                              col_data = mat_all, 
                              overdispersion = fit_rough$overdispersions, 
                              overdispersion_shrinkage = F, 
                              size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                              on_disk = FALSE)
            } else {
                # need a different model formula if only one celltype presents (or split.by = F)
                fit_rough <- glm_gp_disp_only(data = count_data2[idx_for_DE, ], 
                                    design = ~ 1 + log_ct, 
                                    col_data = mat_all, 
                                    size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                    on_disk = FALSE)
                
                fit <- glm_gp(data = count_data2[idx_for_DE, ], 
                              design = ~ 1 + weight + log_ct, 
                              col_data = mat_all, 
                              overdispersion = fit_rough$overdispersions, 
                              overdispersion_shrinkage = F, 
                              size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                              on_disk = FALSE)
            }
            
            # get the SE for each coefficients
            identity_design_matrix <-  diag(nrow = ncol(fit$Beta))
            pred <- predict(fit, se.fit = TRUE, newdata = identity_design_matrix)
            se = pred$se.fit
            
            # get p-val
            beta = fit$Beta
            p = pchisq((beta/se)^2, df = 1, lower.tail = F)
            p = as.data.frame(p)
            beta = as.data.frame(beta)
            names(p) = paste0("p_", names(p))
            names(beta) = paste0("beta_", names(beta))
            p$gene_ID = row.names(p)
            res = cbind(beta, p)
            
            
            ##################################
            ##   2.  LOO DE test 
            for(gene_rm in de_gene_to_rm){
                idx_gene_rm = which(row.names(count_data2) == gene_rm )
                
                if(gene_rm == PRTB){
                    mat_tmp = mat_all[, c("cell_type", "log_ct", "gene")]
                    mat_tmp$weight = 0
                    mat_tmp$weight[mat_tmp$gene != "NT"] = 1
                } else {
                    mat_tmp = mat_all[, c("cell_type", paste0("weight_", gene_rm), "log_ct")]
                    names(mat_tmp) = c("cell_type", "weight", "log_ct")
                }
                
                # ###
                if(length(celltype_list) > 1){
                    fit <- glm_gp(data = count_data2[idx_gene_rm, ], 
                                  design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                                  col_data = mat_tmp, 
                                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                  on_disk = FALSE)
                }else {
                    # need a different model formula if only one celltype presents (or split.by = F)
                    fit <- glm_gp(data = count_data2[idx_gene_rm, ], 
                                  design = ~ 1 + weight + log_ct, 
                                  col_data = mat_tmp, 
                                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                  on_disk = FALSE)
                }
                
                
                # get the SE for each coefficients
                identity_design_matrix <-  diag(nrow = ncol(fit$Beta))
                pred <- predict(fit, se.fit = TRUE, newdata = identity_design_matrix)
                se = pred$se.fit
                
                # get p-val
                beta = fit$Beta
                p = pchisq((beta/se)^2, df = 1, lower.tail = F)
                p = as.data.frame(p)
                beta = as.data.frame(beta)
                names(p) = paste0("p_", names(p))
                names(beta) = paste0("beta_", names(beta))
                p$gene_ID = gene_rm
                tmp_res = cbind(beta, p)
                
                res = rbind(res, tmp_res)
            }
            
            # make idx_for_DE complete again 
            idx_for_DE = c( idx_for_DE, idx_LOO_rm )
            
        } else if (DE_FLAG == "standard"){
            print(paste0(PRTB, " is running standard DE test! "))
            
            # 
            if(length(celltype_list) > 1){
                fit_rough <- glm_gp_disp_only(data = count_data2[idx_for_DE, ], 
                                    design = ~ 0 + cell_type + log_ct, 
                                    col_data = mat_all, 
                                    size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                    on_disk = FALSE)
                
                fit <- glm_gp(data = count_data2[idx_for_DE, ], 
                              design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                              col_data = mat_all, 
                              overdispersion = fit_rough$overdispersions, 
                              overdispersion_shrinkage = F, 
                              size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                              on_disk = FALSE)
            } else {
                # need a different model formula if only one celltype presents (or split.by = F)
                fit_rough <- glm_gp_disp_only(data = count_data2[idx_for_DE, ], 
                                    design = ~ 1 + log_ct, 
                                    col_data = mat_all, 
                                    size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                                    on_disk = FALSE)
                
                fit <- glm_gp(data = count_data2[idx_for_DE, ], 
                              design = ~ 1 + weight + log_ct, 
                              col_data = mat_all, 
                              overdispersion = fit_rough$overdispersions, 
                              overdispersion_shrinkage = F, 
                              size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                              on_disk = FALSE)
            }
            
            # get the SE for each coefficients
            identity_design_matrix <-  diag(nrow = ncol(fit$Beta))
            pred <- predict(fit, se.fit = TRUE, newdata = identity_design_matrix)
            se = pred$se.fit
            
            # get p-val
            beta = fit$Beta
            p = pchisq((beta/se)^2, df = 1, lower.tail = F)
            p = as.data.frame(p)
            beta = as.data.frame(beta)
            names(p) = paste0("p_", names(p))
            names(beta) = paste0("beta_", names(beta))
            p$gene_ID = row.names(p)
            res = cbind(beta, p)
        }
        
        ###########
        #   final step: paste the fold-change info to data.frame "res"
        res = cbind(res, fc_mat[idx_for_DE, ])
        res$DE_method = DE_FLAG
        
        # save the DE results for this PRTB to "all_res"
        all_res[[PRTB]] = res
        
        print(paste0("DE test for '", PRTB, "' is completed. Number of remaining = ", length(all_PRTB_list) - match(PRTB, all_PRTB_list)))
    }
    
    return(all_res)
}



#' Rearrange the DE results into a list of Z-score matrices
#'
#' A function to re-arrange the DE results produced by scoringDE() into a list of Z-score matrices. 
#' Each matrix represents one cell type (if multiple cell types were included), and contains the 
#' DE test Z-scores for each valid gene being tested (rows) and each perturbation target (columns).  
#' 
#' @export
#' 
#' @param de_res the DE results produced by scoringDE(), which is a list of data frames.
#' @param p_threshold the DE P-value threshold to define statistically significant DE genes. 
#' @param fc_threshold the log-fold-change threhsold to define statistically significant DE genes. 
#' 

get_DE_mat = function(de_res = NULL, 
                      p_threshold = 0.05/30000,
                      fc_threshold = 0.2){
    # first get the names of all the perturbations
    PRTB_list = names(de_res)
    
    # and get the names of all the cell types
    celltype_list = grep("fc_", colnames(de_res[[1]]), value = T)
    celltype_list = sort(gsub("fc_", "", celltype_list))
    
    # now we will extract the genes that meet the p_threshold and fc_threshold criteria
    gene_list = c()
    slct_PRTB_list = c()  #  we also filter out PRTB with no responses;
    
    # loop through all the PRTBs
    for(PRTB in PRTB_list){
        # load the results file
        DE_res = de_res[[PRTB]]
        
        isDE = vector()
        for(i in 1:nrow(DE_res)){
            # count if there is any sig DE for each gene being tested
            ct_DE = sum( (DE_res[i, paste0("p_cell_type", celltype_list, ":weight")] <= p_threshold) & 
                             (abs(DE_res[i, paste0("fc_", celltype_list)]) >= fc_threshold), na.rm = T )
            isDE = c(isDE, ct_DE)
        }
        DE_res$ct_DE = isDE
        DE_res = DE_res[order(DE_res$ct_DE, decreasing = T), ]
        
        print(paste(PRTB, ":", sum(DE_res$ct_DE != 0) ) )
        
        if(sum(DE_res$ct_DE != 0) > 0){
            gene_list = c(gene_list, DE_res[DE_res$ct_DE != 0, "gene_ID"])
            slct_PRTB_list = c(slct_PRTB_list, PRTB)
        }
    }
    
    # get the unique gene names
    gene_list2 = unique(gene_list)
    slct_PRTB_list = unique(slct_PRTB_list)
    
    # B. second loop to extract z score based on the DEG list we made in the first loop
    #     here we extract the DEG results for each cell type and make it a DEG matrix individually;
    
    DEG_mat = list() # to save all the DEG matrices 
    logfc_mat = list() # to save the log-fc information
    
    for(CELLTYPE in celltype_list){
        # save the results for each cell type
        merged_dat1= data.frame(gene_ID = gene_list2)
        merged_dat2 = data.frame(gene_ID = gene_list2)
        
        # 
        for(PRTB in slct_PRTB_list){
            DE_res = de_res[[PRTB]]
            
            # 
            DE_res$z = sign(DE_res[, paste0("beta_cell_type", CELLTYPE, ":weight")]) * sqrt(qchisq(DE_res[, paste0("p_cell_type", CELLTYPE, ":weight")], df = 1, lower.tail = F))
            logfc_res = DE_res[, c("gene_ID", paste0("fc_", CELLTYPE))]
            DE_res = DE_res[, c("gene_ID", "z")]
            names(DE_res) = c("gene_ID", PRTB)
            merged_dat1 = merge(merged_dat1, DE_res, by = "gene_ID", all.x = T)
            
            # 
            names(logfc_res) = c("gene_ID", PRTB)
            merged_dat2 = merge(merged_dat2, logfc_res, by = "gene_ID", all.x = T)
            
        }
        # append the matrix to the list 
        DEG_mat[[CELLTYPE]] = merged_dat1
        logfc_mat[[CELLTYPE]] = merged_dat2
        
    }
    
    return(DEG_mat)
}



#' QC for the list of DE Z-score matrices
#' 
#' A function to perform QC/filtering for the DE Z-score matrices that produced by get_DE_mat(). 
#' 
#' @export
#' @param DEG_mat the list of DE Z-score matrices that produced by get_DE_mat()
#' @param zscore_cap the cap value for the Z-scores. Any absolute(Z-score) larger than this value will be 
#' set to this value to avoid extreme values affecting the downstream analyses. Default is 37 
#' which is the machine precision limit for a Z-score to produce a non-zero P-value (= ~1e-300). 
#' @param mask_target a boolen value to indicate if the Z-score of the perturbation target (labelled 
#' by the column names) should be masked as 0. Default is FALSE so no masking will happen.
#' @param min_sig_DEG the minimal number of significant DE genes that each column must contain. Any 
#' column with sig DE genes less than this value will be removed from the matrix. Default is 0 so no 
#' column will be removed. 
#' @param p_threshold the P-value threshold to define the 'significant' DE genes. Default is 0.05/30000, which
#' is approximately the Bonferroni correction threshold for genome-wide DE tests (assuming 30,000 genes in total).
#' @param center a boolen value to indicate whether we should center each column to 0. Default is FALSE.
#' @param scale a boolen value to indicate whether we should scale each column to have variance = 1. Default is FALSE.
#' 
#' @return The function will return a list of QCed DE Z-score matrices. This list can directly be the input for 
#' DEmultiCCA(). 
#' 

prune_DE_mat = function(DEG_mat = NULL, 
                        zscore_cap = 37, 
                        mask_target = FALSE, 
                        p_threshold = 0.05/30000,
                        min_sig_DEG = 0, 
                        center = FALSE, 
                        scale = FALSE,
                        ...){
    # get the cell type names and the column names (list of PRTBs) from the DEG_mat 
    celltype_list = names(DEG_mat)
    PRTB_list = colnames(DEG_mat[[1]])[-1]
    
    # get the names for all the DE genes
    gene_ID = DEG_mat[[1]]$gene_ID
    
    # If TRUE, we will mask the PRTB target in each PRTB to 0:
    if(mask_target == TRUE){
        for(CELLTYPE in celltype_list){
            tmp = DEG_mat[[CELLTYPE]]
            for(PRTB in PRTB_list){
                if(PRTB %in% tmp$gene_ID){
                    tmp[tmp$gene_ID %in% PRTB, PRTB] = 0
                } 
            }
            DEG_mat[[CELLTYPE]] = tmp
        }
    }

    # reset all the NAs to 0s
    for(i in 1:length(DEG_mat)){
        DEG_mat[[i]][is.na(DEG_mat[[i]])] = 0
    }
    
    # we also cap the z-score if a zscore_cap is provided. 
    for(i in 1:length(DEG_mat)){
        DEG_mat[[i]][, -1][DEG_mat[[i]][, -1] >= zscore_cap] = zscore_cap
        DEG_mat[[i]][, -1][DEG_mat[[i]][, -1] <= -zscore_cap] = -zscore_cap
    }
    
    # remove any column with all 0s (this will cause an error when doing SVD)
    DEG_mat2 = list()
    for(CELLTYPE in celltype_list){
        tmp = DEG_mat[[CELLTYPE]][, -1, drop = F]
        rm_idx = vector()
        for(i in 1:ncol(tmp)){
            if(sum(abs(tmp[, i]) > sqrt(qchisq(p_threshold, df=1, lower.tail = F))) < min_sig_DEG){
                rm_idx = c(rm_idx, i)
            }
        }
        if(length(rm_idx) != 0){
            tmp = tmp[, -rm_idx, drop = F]
        }
        row.names(tmp) = DEG_mat[[CELLTYPE]][, 1]
        # colnames(tmp) = tmp[, "gene_ID"]
        if(ncol(tmp) > 1){
            DEG_mat2[[CELLTYPE]] = tmp
        } else {
            DEG_mat2[[CELLTYPE]] = NULL
        }
    }
    
    # centering and scaling 
    if(center == TRUE){
        for(i in 1:length(DEG_mat2)){
            DEG_mat2[[i]] = scale(DEG_mat2[[i]], center = center, scale = scale)
        }
    }
    # 
    return(DEG_mat2)
}












