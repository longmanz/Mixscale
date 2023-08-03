#' @include get_fold_change.R
#'
NULL

#' Scoring-based DE test
#'
#' Function to perform differential expression (DE) tests based on the perturbation scores from the 
#' PRTBScoring() function. It is a multivariate negative binomial based model that incorporates both the heterogeneity 
#' of perturbation strength in each cell, as well as their cell type background. 
#'
#' @inheritParams Seurat::FindMarkers
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene labels.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
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


scoringDE = function (object, assay = "PRTB", slot = "data", labels = "gene", 
                      nt.class.name = "NT", verbose = FALSE, 
                      split.by = "cell_type",
                      total_ct_labels = "nCount_RNA",  
                      logfc.threshold = 0.2, 
                      ...) {
    # 
    print("Running scoring DE test!\n")
    
    # check if the required package is installed
    glmGamPoi.installed <- PackageCheck("glmGamPoi", error = FALSE)
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
    
    dat <- GetAssayData(object = object[[assay]], 
                        slot = "data")[de.genes, all.cells, drop = FALSE]
    if (slot == "scale.data") {
        dat <- ScaleData(object = dat, features = de.genes, 
                         verbose = FALSE)
    }
    
    # get the PRTBs with weights
    DefaultAssay(object) = "PRTB"
    prtb_score <- Tool(object = object, slot = "PRTBScoring")
    wt_PRTB_list = sort(names(prtb_score))
    DefaultAssay(object) = "RNA"
    
    # get the full list of PRTBs 
    all_PRTB_list = sort(unique(object[[labels]]))
    all_PRTB_list = all_PRTB_list[all_PRTB_list != nt.class.name]
    wt_PRTB_list = wt_PRTB_list[wt_PRTB_list %in% all_PRTB_list]
    
    # 
    # count_data = GetAssayData(object = object[['RNA']], slot = "counts")   # get the raw count data (required by neg binom)
    # count_data_std = GetAssayData(object = object[['RNA']], slot = "data")   # get the std data (required by log-fold change test)
    
    # get the overall covariate matrix from the meta.data
    mat_B = data.frame(cell_label = colnames(object), 
                       nCount_RNA = object[[total_ct_labels]], 
                       cell_type = as.factor(object[[split.by]]), 
                       gene = object[[labels]] )
    
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

            # DE Flag:
            DE_FLAG = "standard"
        }
        
        ##############################################
        # get idx for the cells of this PRTB
        idx = match(mat_all$cell_label, colnames(object))

        idx_NT = match(mat_all$cell_label[mat_all$gene == "NT"], colnames(count_data2))
        idx_P = match(mat_all$cell_label[mat_all$gene != "NT"], colnames(count_data2))

        count_data2 = GetAssayData(object = object[['RNA']], slot = "counts")[, idx]
        count_data_std2 = GetAssayData(object = object[['RNA']], slot = "data")[, idx]

        # do not do fold-change check, only do var and min.pct check
        idx_for_DE = which(apply(count_data_std2, MARGIN = 1, FUN = get_idx, idx_P = idx_P, idx_NT = idx_NT, logfc.threshold = 0))
        # # always include the PRTB gene itself in the DE test
        # idx_PRTB = which(row.names(count_data_std2) %in% PRTB)
        # # merge the two indices
        # idx_for_DE = sort(unique(c(idx_for_DE, idx_PRTB)))
        
        # here we will calculate the logfc within each cell type. 
        for( idx_i in 1:(length(celltype_list)) ){
            celltype = celltype_name = levels(mat_all$cell_type)[idx_i]
            
            # select the count matrix for one celltype 
            idx = match(mat_all[ mat_all$celltype == celltype, "cell_label"], colnames(count_data2))
            count_data_celltype = count_data2[idx_for_DE, idx] # make sure to include idx_for_DE
            count_data_std_celltype = count_data_std2[idx_for_DE, idx]
            
            # 
            idx_NT_celltype = match(mat_all$cell_label[mat_all$celltype == celltype & 
                                                           mat_all$gene == "NT"], colnames(count_data_celltype))
            idx_P_celltype = match(mat_all$cell_label[mat_all$celltype == celltype & 
                                                          mat_all$gene != "NT"], colnames(count_data_celltype))
            # do not do fold-change check, only do var and min.pct check
            fc = apply(count_data_std_celltype, MARGIN = 1, FUN = get_fc, 
                       idx_P = idx_P_celltype, idx_NT = idx_NT_celltype)
            # idx_P = c(idx_P_celltype, idx_NP_celltype), idx_NT = idx_NT_celltype)
            
            res[paste0("fc_", celltype_name)] =  fc
            rm(fc)
        }
        
        
        
    }
    
}
    






