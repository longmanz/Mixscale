#'
NULL


#' Mixscale scoring for perturbations
#'
#' Function to calculate perturbation scores for perturbed and non-perturbed gRNA expressing cells.
#' The perturbation score reflects the perturbation strength of each cells (inherited from the RunMixscape()
#' function). It is calculated by using the large-effect DE genes from raw DE tests between the 
#' perturbed and non-perturbed gRNA expressing cells. 
#'
#' @export
#' 
#' @inheritParams Seurat::RunMixscape
#' @import Seurat
#' 
#' @param object An object of class Seurat.
#' @param assay Assay to use for mixscape classification.
#' @param slot Assay data slot to use.
#' @param labels metadata column with target gene labels.
#' @param nt.class.name Classification name of non-targeting gRNA cells.
#' @param new.class.name Name of mixscale scores to be stored in
#' metadata.
#' @param min.de.genes Required number of genes that are differentially
#' expressed for method to separate perturbed and non-perturbed cells.
#' @param min.cells Minimum number of cells in target gene class. If fewer than
#' this many cells are assigned to a target gene class during classification,
#' all are assigned NP.
#' @param de.assay Assay to use when performing differential expression analysis.
#' Usually RNA.
#' @param logfc.threshold the log-fold-change threshold to select the large-effect
#' DE genes. Only DE genes with log-fold-change larger than this value will be 
#' selected. Default is 0.25.
#' @param verbose Display messages
#' @param split.by metadata column with experimental condition/cell type
#' classification information. This is meant to be used to account for cases a
#' perturbation is condition/cell type -specific.
#' @param fine.mode When this is equal to TRUE, DE genes for each target gene
#' class will be calculated for each gRNA separately and pooled into one DE list
#' for calculating the perturbation score of every cell and their subsequent
#' classification.
#' @param fine.mode.labels metadata column with gRNA ID labels.
#' @param DE.gene specify a list of user-defined large-effect DE genes to calculate the perturbation score.
#' @param max.de.genes the maximum number of top large-effect DE genes to calculate the perturbation score. Default is 100. 
#' @param harmonize a boolen value to specify whether a harmonization of the cell-type proportion between the NT cells and 
#' the perturbed cells should be performed prior to the DE test. If fine.mode is TRUE, this harmonization step will be 
#' performed for each fine.mode gRNA. Default is FALSE. 
#' @param min_prop_ntgd a minimal threshold to remove cells if any cell type has a proportion less than this value. It will
#' only be used when harmonize is TRUE. Default is 0.1. 
#' @param pval.cutoff specify the DE test p-value cutoff (after Bonferroni correction) to select top large-effect DE genes.
#' Default is 0.05. 
#' 
#' @return Returns a Seurat object containing the perturbation scores. It is stored in the Tool Data of the object, not in the 
#' meta.data. 
#' @concept perturbation_scoring

RunMixscale = function (object, assay = "PRTB", slot = "scale.data", labels = "gene", 
                        nt.class.name = "NT", 
                        new.class.name = "mixscale_score", 
                        min.de.genes = 5, min.cells = 5, de.assay = "RNA", logfc.threshold = 0.25, 
                        verbose = FALSE, split.by = NULL, fine.mode = FALSE, 
                        fine.mode.labels = "guide_ID", 
                        DE.gene = NULL, max.de.genes = 100, harmonize = F, 
                        min_prop_ntgd = 0.1, pval.cutoff = 0.05, 
                        seed = 10282021) 
{
    print("Running Mixscale scoring to calculate the perturbation scores \n")
    
    assay <- assay %||% DefaultAssay(object = object)
    if(!assay %in% names(object@assays)){
        stop(paste0("The 'assay' being specified does not exist! Please check. Have you run CalcPerturbSig() yet?"))
    }
    
    if (is.null(x = labels)) {
        stop("Please specify target gene class metadata name")
    }
    prtb_markers <- list()
    prtb_markers2 <- list()
    object[[new.class.name]] <- object[[labels]]
    object[[new.class.name]][, 1] <- as.character(x = object[[new.class.name]][, 
                                                                               1])
    gv.list <- list()
    if (is.null(x = split.by)) {
        split.by <- splits <- "con1"
    } else {
        splits <- as.character(x = unique(x = object[[split.by]][, 
                                                                 1]))
    }
    cells.s.list <- list()
    
    # 
    Idents(object = object) <- "con1"
    cells.s <- WhichCells(object = object, idents = "con1")
    # cells.s.list[[s]] <- cells.s
    genes <- setdiff(x = unique(x = object[[labels]][cells.s, 
                                                     1]), y = nt.class.name)
    # 
    for (gene in genes) {
        Idents(object = object) <- labels
        
        if (isTRUE(x = verbose)) {
            message("Processing ", gene)
        }
        orig.guide.cells <- intersect(x = WhichCells(object = object, 
                                                     idents = gene), y = cells.s)
        nt.cells <- intersect(x = WhichCells(object = object, 
                                             idents = nt.class.name), y = cells.s)
        
        #############################################
        
        if (isTRUE(x = fine.mode)) {
            guides <- setdiff(x = unique(x = object[[fine.mode.labels]][orig.guide.cells, 
                                                                        1]), y = nt.class.name)
            all.de.genes <- c()
            for (gd in guides) {
                gd.cells <- rownames(x = object[[]][orig.guide.cells, 
                ])[which(x = object[[]][orig.guide.cells, 
                                        fine.mode.labels] == gd)]
                # we will need to extract the NT cells based on each celltype and do harmonization based on cell comp in PRTB cells
                if(harmonize == T){
                    # this is a flag to indicate if the harmonization process is okay (selected_Cell >= 50% of total cell)
                    flag_good_harm = F
                    # the initial list of split-groups
                    splits_list = splits 
                    
                    while(flag_good_harm == F){
                        # get the cell label splitted by splits
                        cells.s.list.gd = list()
                        cells.s.list.ntgd = list()
                        Idents(object = object) <- split.by
                        
                        for (s in splits_list) {
                            cells.s.list.gd[[s]] <- intersect(gd.cells, WhichCells(object = object, idents = s))
                            cells.s.list.ntgd[[s]] <- intersect(nt.cells, WhichCells(object = object, idents = s))
                        }
                        
                        # calculate the desired number of nt cells in each splits;
                        length.gd = sapply(X = cells.s.list.gd, FUN = length)
                        length.ntgd = sapply(X = cells.s.list.ntgd, FUN = length)
                        
                        prop.gd = length.gd/sum(length.gd, na.rm = T)
                        # prop.ntgd = length.ntgd/sum(length.ntgd, na.rm = T)
                        sum.desire.length.ntgd = floor(min(length.ntgd/prop.gd, na.rm = T))
                        if(sum.desire.length.ntgd > sum(length.ntgd, na.rm = T)){
                            stop("The sum.desire.length.ntgd is greater than the total number of NT cells. Need to check!")
                        }
                        
                        if(sum.desire.length.ntgd >= min_prop_ntgd*sum(length.ntgd, na.rm = T) ){
                            flag_good_harm = T
                        } else {
                            print(paste("Removing cell from ", splits_list[which.min(length.ntgd)], "due to 50% check during harmonization step."))
                            splits_list = splits_list[-which.min(length.ntgd)]
                        }
                    }
                    
                    # calculate the final number of NT cells to extract from splits_list
                    desire.length.ntgd = floor(sum.desire.length.ntgd*prop.gd)
                    
                    # start to subsample the nt cells based on the desire length:
                    sub.cells.s.list.ntgd = list()
                    for (s in splits_list) {
                        set.seed(seed = seed)
                        sub.cells.s.list.ntgd[[s]] <- sample(x = cells.s.list.ntgd[[s]], size = desire.length.ntgd[s])
                    }
                    
                    # collapse the list into a single vectors of sub-sampled NT cells
                    sub.ntgd.cells = Reduce(c, sub.cells.s.list.ntgd)
                    rm(cells.s.list.gd, cells.s.list.ntgd, length.gd, length.ntgd, prop.gd, sum.desire.length.ntgd, desire.length.ntgd, sub.cells.s.list.ntgd)
                    if(verbose){
                        print("Done with harmonizing the cell composition in NT cells (fine mode).")
                    }
                    
                    Idents(object = object) <- labels
                } else {
                    sub.ntgd.cells = nt.cells
                }
                
                # run DE
                if(!is.null(DE.gene) ){
                    if(!is.null(DE.gene[[gene]]) | length(DE.gene[[gene]]) != 0){
                        all.de.genes = DE.gene[[gene]]
                    } else {
                        all.de.genes = character()
                        print(paste("No de.genes are provided for PRTB:", gene, ". Pls check!"))
                    }
                } else {
                    # run DE 
                    de.genes <- Seurat:::TopDEGenesMixscape(object = object, 
                                                   ident.1 = gd.cells, ident.2 = sub.ntgd.cells, de.assay = de.assay, 
                                                   logfc.threshold = logfc.threshold, labels = fine.mode.labels, 
                                                   verbose = verbose, pval.cutoff = pval.cutoff)
                    all.de.genes <- c(all.de.genes, de.genes)
                }
            }
            all.de.genes <- unique(all.de.genes)
        } else {
            # we will need to extract the NT cells based on each celltype and do harmonization based on cell comp in PRTB cells
            if(harmonize == T){
                # this is a flag to indicate if the harmonization process is okay (selected_Cell >= 50% of total cell)
                flag_good_harm = F
                # the initial list of split-groups
                splits_list = splits 
                
                while(flag_good_harm == F){
                    # get the cell label splitted by splits
                    cells.s.list.gene = list()
                    cells.s.list.nt = list()
                    Idents(object = object) <- split.by
                    
                    for (s in splits_list) {
                        cells.s.list.gene[[s]] <- intersect(orig.guide.cells, WhichCells(object = object, idents = s))
                        cells.s.list.nt[[s]] <- intersect(nt.cells, WhichCells(object = object, idents = s))
                    }
                    
                    # calculate the desired number of nt cells in each splits;
                    length.gene = sapply(X = cells.s.list.gene, FUN = length)
                    length.nt = sapply(X = cells.s.list.nt, FUN = length)
                    
                    prop.gene = length.gene/sum(length.gene, na.rm = T)
                    # prop.nt = length.nt/sum(length.nt, na.rm = T)
                    sum.desire.length.nt = floor(min(length.nt/prop.gene, na.rm = T))
                    if(sum.desire.length.nt > sum(length.nt, na.rm = T)){
                        stop("The sum.desire.length.nt is greater than the total number of NT cells. Need to check!")
                    }
                    # 
                    if(sum.desire.length.nt >= min_prop_ntgd*sum(length.nt, na.rm = T) ){
                        flag_good_harm = T
                    } else {
                        print(paste("Removing cell from ", splits_list[which.min(length.nt)], "due to 50% check during harmonization step."))
                        splits_list = splits_list[-which.min(length.nt)]
                    }
                }
                
                ###
                desire.length.nt = floor(sum.desire.length.nt*prop.gene)
                
                # start to subsample the nt cells based on the desire length:
                sub.cells.s.list.nt = list()
                for (s in splits_list) {
                    set.seed(seed = seed)
                    sub.cells.s.list.nt[[s]] <- sample(x = cells.s.list.nt[[s]], size = desire.length.nt[s])
                }
                
                # collapse the list into a single vectors of sub-sampled NT cells
                sub.nt.cells = Reduce(c, sub.cells.s.list.nt)
                rm(cells.s.list.gene, cells.s.list.nt, length.gene, length.nt, prop.gene, sum.desire.length.nt, desire.length.nt, sub.cells.s.list.nt)
                if(verbose){
                    print("Done with harmonizing the cell composition in NT cells.")
                }
                
                Idents(object = object) <- labels
            } else {
                sub.nt.cells = nt.cells
            }
            # run DE
            if(!is.null(DE.gene) ){
                if(!is.null(DE.gene[[gene]]) | length(DE.gene[[gene]]) != 0){
                    all.de.genes = DE.gene[[gene]]
                } else {
                    all.de.genes = character()
                    print(paste("No de.genes are provided for PRTB:", gene, ". Pls check!"))
                }
            } else {
                all.de.genes <- Seurat:::TopDEGenesMixscape(object = object, 
                                                   ident.1 = orig.guide.cells, ident.2 = sub.nt.cells, 
                                                   de.assay = de.assay, logfc.threshold = logfc.threshold, 
                                                   labels = labels, verbose = verbose, pval.cutoff = pval.cutoff)
            }
            
            
        }
        # print(gene)
        # print(all.de.genes)
        
        # only keep the top max.de.genes as the de.genes for PRTB score calculation
        if(!is.null(max.de.genes) ){
            if(length(all.de.genes) <= max.de.genes){
                if(verbose){
                    message(paste("The number of de.genes (", length(all.de.genes), ") is less than max.de.genes (", max.de.genes, ")."))
                }
            } else {
                if(verbose){
                    message(paste("The number of de.genes (", length(all.de.genes), ") is larger than max.de.genes (", max.de.genes, ").",  
                                "Restricting to top", max.de.genes, "genes..."))
                }
                all.de.genes = all.de.genes[1:max.de.genes]
            }
        }
        
        # use user-defined DE.gene list
        # if(!is.null(DE.gene) ){
        #     all.de.genes = DE.gene
        # }
        
        for (s in splits) {
            Idents(object = object) <- split.by
            cells.s.list[[s]] <- WhichCells(object = object, idents = s)
            
            prtb_markers[[s]][[gene]] <- all.de.genes
            prtb_markers2[[s]][[gene]] <- all.de.genes
            if (length(x = all.de.genes) < min.de.genes) {
                prtb_markers[[s]][[gene]] <- character()
            }
            
        }
        if(verbose){
            message(paste0("Done with extracting top DE genes for ", gene))
        }
    }
    
    all_markers <- unique(x = unlist(x = prtb_markers))
    missing_genes <- all_markers[!all_markers %in% rownames(x = object[[assay]])]
    # print(missing_genes)
    object <- Seurat:::GetMissingPerturb(object = object, assay = assay, 
                                         features = missing_genes, verbose = verbose)
    
    if(verbose){
        message("Done with getting Missing PRTB")
    }
    
    for (s in splits) {
        print(splits)
        cells.s <- cells.s.list[[s]]
        genes <- setdiff(x = unique(x = object[[labels]][cells.s, 
                                                         1]), y = nt.class.name)
        for (gene in genes) {
            Idents(object = object) <- labels
            post.prob <- 0
            orig.guide.cells <- intersect(x = WhichCells(object = object, 
                                                         idents = gene), y = cells.s)
            nt.cells <- intersect(x = WhichCells(object = object, 
                                                 idents = nt.class.name), y = cells.s)
            all.cells <- c(orig.guide.cells, nt.cells)
            if (length(x = prtb_markers[[s]][[gene]]) == 0) {
                if (verbose) {
                    message("  Fewer than ", min.de.genes, " DE genes for ", 
                            gene, ". Assigning cells as NP.")
                }
                object[[new.class.name]][orig.guide.cells, 1] <- paste0(gene, " NP")
            }
            else {
                if (verbose) {
                    message("  ", gene)
                }
                de.genes <- prtb_markers[[s]][[gene]]
                dat <- GetAssayData(object = object[[assay]], 
                                    slot = "data")[de.genes, all.cells, drop = FALSE]
                if (slot == "scale.data") {
                    dat <- ScaleData(object = dat, features = de.genes, 
                                     verbose = FALSE)
                }

                # the first step to calculate the overall PRTB score
                if(verbose){
                    cat(paste0("Calculating the overall PRTB score...\n"))
                }
                Idents(object = object) <- new.class.name
                guide.cells <- intersect(x = WhichCells(object = object, 
                                                        idents = gene), y = cells.s)
                vec <- matrixStats::rowMeans2(x = dat[, guide.cells, drop = FALSE]) - 
                    matrixStats::rowMeans2(x = dat[, nt.cells, drop = FALSE])
                
                # save the mat and vec to easily calculate the weights by rowSums
                pvec_mat = sweep(t(dat), MARGIN=2, vec, `*`) 
                vec_mat = vec * vec
                names(vec_mat) = colnames(pvec_mat)
                
                # the weights 
                pvec = matrixStats::rowSums2(pvec_mat)/sum(vec_mat)
                names(pvec) = rownames(pvec_mat)
                
                # create a list to store the PRTB score
                gv <- as.data.frame(x = pvec)
                gv[, labels] <- nt.class.name
                gv[intersect(x = rownames(x = gv), y = guide.cells), 
                   labels] <- gene
                gv.list[[gene]][[s]] <- gv
                
                # the LOOv2 weights
                # for(omit_gene in de.genes){
                #     if(verbose){
                #         cat(paste0("Calculating the LOO PRTB score by using rowSums2 for ", omit_gene, " in subset ", s, "...\n"))
                #     }
                #     remain_gene = de.genes[which(de.genes != omit_gene)]
                #     pvec2 <- rowSums2(pvec_mat[, remain_gene, drop = F])/sum(vec_mat[remain_gene])
                #     # save the LOO weights 
                #     gv.list[[gene]][[s]][, omit_gene] <- pvec2
                # }
                
                # 2023 July 03: substitute the above loop with the following matrix manipulation for speed. 
                omit_mat <- outer(de.genes, de.genes, `!=`)
                
                # Function to calculate pvec2
                calc_pvec2 <- function(include_gene) {
                    pvec2 <- matrixStats::rowSums2(pvec_mat[, include_gene, drop = F])/sum(vec_mat[include_gene])
                    return(pvec2)
                }
                
                # Apply the function to each column of omit_mat
                gv.list[[gene]][[s]][, de.genes] <- apply(omit_mat, 2, calc_pvec2)
                
                
                # the second step to calculate the leave-one-out (LOO) PRTB score  
                if(verbose){
                    cat(paste0("Done calculating LOO PRTB score for ", length(de.genes), " genes in ", s, "...\n"))
                }

            }
            
        }
        message(paste0("Done with calculating scores for ", s))
    }
    # SeuratObject::Tool(object = object) <- gv.list
    
    # added Jan 16: calculate the standardized scores and append them to the meta-data
    # get the list of PRTBs 
    wt_PRTB_list = sort(names(gv.list))
    all_PRTB_list = sort(unique(object[[labels]][,1]))
    all_PRTB_list = all_PRTB_list[all_PRTB_list != nt.class.name]
    wt_PRTB_list = wt_PRTB_list[wt_PRTB_list %in% all_PRTB_list]
    
    # 
    if(is.null(split.by)) {
        mat_B = data.frame(cell_label = colnames(object), 
                           gene = object[[labels]][,1] )
    } else {
        mat_B = data.frame(cell_label = colnames(object), 
                           gene = object[[labels]][,1] )
    }
    
    # 
    all_score = data.frame() # to store the scores from each PRTB
    for(PRTB in all_PRTB_list){
        mat_A = data.frame()
        # check if the scores are calculated successful for this PRTB
        if(PRTB %in% wt_PRTB_list){
            celltype_list = names(gv.list[[PRTB]])
            for(celltype in celltype_list){
                tmp = gv.list[[PRTB]][[celltype]][, c("pvec", "gene"), drop = FALSE]

                # get the idx for NT cells and PRTBed cells
                idx_NT = which(tmp$gene == nt.class.name)
                idx_gene = which(tmp$gene == PRTB)
                
                #  1. calculate the overall weights
                # calculate the mean and sd of the PRTB score for the NT cells
                mean_NT = mean(tmp$pvec[idx_NT], na.rm = T)
                sd_NT = sd(tmp$pvec[idx_NT], na.rm = T)
                # standardize the PRTB scores for the PRTBed cells based on the mean and SD from those of the NT cells
                std_weight_gene = (tmp$pvec[idx_gene] - mean_NT)/sd_NT
                # convert those negative standardised PRTB score to 0
                # std_weight_gene[which(std_weight_gene < 0)] = 0
                
                # create a new column called "weight" in the tmp dataframe. 
                tmp$weight = 0 
                tmp$weight[idx_gene] = std_weight_gene
                
                tmp$cell_label = row.names(tmp)
                tmp = tmp[, c("cell_label", "gene", "pvec", "weight")]
                
                mat_A = rbind(mat_A, tmp)
                rm(tmp)
            }
        } else {
            celltype_list = names(gv.list[[1]])
            # 
            tmp = mat_B[mat_B$gene %in% c(PRTB, nt.class.name), ]
            tmp$weight = 0
            tmp[tmp$gene == PRTB, "weight"] = 1
            tmp$pvec = tmp$weight
            # 
            tmp = tmp[, c("cell_label", "gene", "pvec", "weight")]
            mat_A = tmp
            rm(tmp)
        }
        all_score = rbind(all_score, mat_A)
    }
    
    # some final editing
    all_score = all_score[!duplicated(all_score$cell_label), ]
    rownames(all_score) = all_score$cell_label
    all_score = all_score[, "weight", drop = FALSE]
    names(all_score) = new.class.name
    
    # add the standardized scores to the meta-data
    object = AddMetaData(object, metadata = all_score)
    
    return(object)
}

