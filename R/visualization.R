#' @include enrichment_test.R
#' @include decomposition.R
NULL


# This script contains multiple functions for visualization of different output in 
# this package. The major part is for the visualization for the decomposition.R. The 
# other part also contains the visualization for the enrichment_test.R.


#' Plot to show the distribution of perturbation scores 
#' 
#' This function will generate a density (ridge) plot for the perturbation scores across
#' different cell types or different perturbation targets. 
#' 
#' @export
#' @import ggplot2
#' @import Seurat
#' @import ggridges
#' 
#' @param object a seurat object returned by RunMixscale()
#' @param labels the column name in the object's meta.data that contains the target
#' gene labels
#' @param nt.class.name the classification name of non-targeting gRNA cells
#' @param split.by metadata column with experimental condition/cell type classification 
#' information. This is meant to be used to account for cases a perturbation is 
#' condition/cell type -specific.
#' @param PRTB the perturbation target genes to extract for plotting. Multiple values are 
#' allowed.
#' @param slct_split.by if only a subset of the conditions/cell-types in the split.by column need 
#' to be plotted, users can specify them as a character vector here. Default is NULL, meaning all the 
#' conditions/cell-types need to be plotted.
#' @param facet_wrap whether to divide the plot into multiple facets based on either the 
#' perturbation targets ("gene") or conditions/cell types ("split.by"). Default is NULL, meaning
#' no facet.
#' @param facet_scale whether to use a fixed scale for y-axis across all facets or allow 
#' y axis to vary. 
#' @param facet_nrow the number of rows to plot the different panels when facet_wrap is set. 
#' 
#' 
#' @return a ggplot2 object that contains the ridge plot.

Mixscale_RidgePlot = function(object = NULL, 
                               labels = "gene", 
                               nt.class.name = "NT", 
                               split.by = NULL, 
                               PRTB = NULL,
                               slct_split.by = NULL,
                               facet_wrap = c(NULL, "gene", "split.by"),
                               facet_scale = c("fixed", "free_y"),
                               facet_nrow = 1, 
                               ...){
    # first get the full object of the PRTB scores
    prtb_score <- Tool(object = object, slot = "RunMixscale")
    # check if the scores exist
    if(is.null(prtb_score)){
        stop(paste0("It seems the scores are not calculated for this object. Please check!"))
    }
    
    # if the split.by is set, need to extract it 
    if (is.null(x = split.by)) {
        split.by <- splits <- "con1"
    } else {
        splits <- as.character(x = unique(x = object[[split.by]][, 
                                                                 1]))
    }
    
    # if slct_split.by is set, we restrict it only to those ones
    if(!is.null(x = slct_split.by)){
        splits = intersect(slct_split.by, splits)
        if(length(splits) == 0){
            stop("The selected slct_split.by has no intersection with the split.by column. Please check.")
        }
    }
    
    # to extract useful info for each PRTB 
    all_scores = data.frame()
    # 
    for (CELLTYPE in splits){
        for (prtb in PRTB){
            if(! prtb %in% names(prtb_score)){
                print(paste0("The provided perturbation: ", prtb, " does not have perturbation scores. Skipping..."))
                next()
            }
            
            scores = prtb_score[[prtb]][[CELLTYPE]][, c(1:2)]  # get the column 1 (score) and 2 (label)
            scores$cell_ID = rownames(scores)
            scores$celltype = CELLTYPE
            scores$PRTB_group = prtb
            
            all_scores = rbind(all_scores, scores)
        }
    }
    
    # adding some more columns 
    all_scores$status = nt.class.name
    all_scores$status[all_scores[[labels]] != nt.class.name] = "perturbed"
    
    all_scores$celltype = factor(x = all_scores$celltype, levels = splits)
    all_scores$status = factor(all_scores$status, levels = c(nt.class.name, "perturbed"))
    
    # generate the ggplot2 object
    if(facet_wrap == "split.by"){
        p3 = ggplot(all_scores, aes(x = pvec, y = celltype, fill = status)) +
            ggridges::geom_density_ridges(scale = 1.4, rel_min_height = 0.01) +
            ggridges::theme_ridges() + xlab("PRTB score") + ylab("condition") +
            facet_wrap(~ PRTB_group, nrow = facet_nrow, scales = facet_scale)
    } else if(facet_wrap == "gene"){
        p3 = ggplot(all_scores, aes(x = pvec, y = PRTB_group, fill = status)) +
            ggridges::geom_density_ridges(scale = 1.4, rel_min_height = 0.01) +
            ggridges::theme_ridges() + xlab("PRTB score") + ylab("gene") +
            facet_wrap(~ celltype, nrow = facet_nrow, scales = facet_scale)
    } else {
        p3 = ggplot(all_scores, aes(x = pvec, y = PRTB_group, fill = status)) +
            ggridges::geom_density_ridges(scale = 1.4, rel_min_height = 0.01) +
            ggridges::theme_ridges() + xlab("PRTB score") + ylab("gene") 
    }
    return(p3)
}





#' Plot to compare the expression level of the perturbation target and the 
#' perturbation scores
#' 
#' This function will generate a connected scatterplot to compare the mean
#' expression level of the perturbation target gene within different perturbation 
#' percentile bins. After running the RunMixscale() function, user can specify the 
#' gene name of the perturbation target and the number of bins to divide the scores 
#' into, and this function will sutomatically generate a connected scatterplot. 
#' Multiple perturbation targets and cell types are allowed. 
#' 
#' @export
#' @import ggplot2
#' @import Seurat
#' 
#' @param object a seurat object returned by RunMixscale()
#' @param assay the assay name to extract the expression level data from for plotting
#' @param slot the slot name to extract the expression level data from for plotting
#' @param labels the column name in the object's meta.data that contains the target
#' gene labels
#' @param nt.class.name the classification name of non-targeting gRNA cells
#' @param split.by metadata column with experimental condition/cell type classification 
#' information. This is meant to be used to account for cases a perturbation is 
#' condition/cell type -specific.
#' @param PRTB the perturbation target genes to extract for plotting. Multiple values are 
#' allowed.
#' @param nbin the number of bins to divide the perturbation scores into. 
#' @param facet_wrap whether to divide the plot into multiple facets based on either the 
#' perturbation targets ("gene") or conditions/cell types ("split.by"). Default is NULL, meaning
#' no facet.
#' @param facet_scale whether to use a fixed scale for y-axis across all facets or allow 
#' y axis to vary. 
#' 
#' @return a ggplot2 object that contains the connected scatterplot. 
#' 
#' 

Mixscale_ScatterPlot = function(object = NULL, 
                                 assay = "RNA", 
                                 slot = "data", 
                                 labels = "gene", 
                                 nt.class.name = "NT", 
                                 split.by = NULL, 
                                 PRTB = NULL, 
                                 nbin = 10,
                                 facet_wrap = c(NULL, "gene", "split.by"),
                                 facet_scale = c("fixed", "free_y"),
                                 ...){
    # if the split.by is set, need to extract it 
    if (is.null(x = split.by)) {
        split.by <- splits <- "con1"
    } else {
        splits <- as.character(x = unique(x = object[[split.by]][, 
                                                                 1]))
    }
    
    # to extract useful info for each PRTB 
    final_res = data.frame()
    # get the full object of the PRTB scores
    prtb_score <- Tool(object = object, slot = "RunMixscale")
    # check if the scores exist
    if(is.null(prtb_score)){
        stop(paste0("It seems the scores are not calculated for this object. Please check!"))
    }
    
    # 
    for (CELLTYPE in splits){
        for (prtb in PRTB){
            if(! prtb %in% names(prtb_score)){
                print(paste0("The provided perturbation: ", prtb, " does not have perturbation scores. Skipping..."))
                next()
            }
            
            scores = prtb_score[[prtb]][[CELLTYPE]][, c(1:2)]  # get the column 1 (score) and 2 (label)
            scores$cell_ID = rownames(scores)
            
            #  get the target gene expression level
            target_expression <- GetAssayData(object = object, assay = assay, slot = slot)[prtb, scores$cell_ID]
            
            target_expression = as.data.frame(target_expression)
            target_expression$cell_ID = rownames(target_expression)
            
            # merge the target expression with prtb score
            combine_dat = merge(scores, target_expression, by = "cell_ID", all.x = T)
            
            # we will keep a record of the average log(expression-level) of the target in the NT cells 
            mean_log_exp_NT = mean(combine_dat[combine_dat[[labels]] == nt.class.name, "target_expression"])
            
            # only focus on PRTB cells and remove NT cell
            combine_dat = combine_dat[combine_dat[[labels]] != nt.class.name, ]
            
            # calculate the cutoff for each bin of the scores
            bin_cutoff = quantile(combine_dat$pvec, probs = seq(0, 1, 1/nbin))
            
            # we will store the expression level of the target gene in each pvec bin
            exp_each_bin = vector()
            
            for(i in 1:nbin){
                cutoff1 = bin_cutoff[i]
                cutoff2 = bin_cutoff[i + 1]
                # 
                idx = which(combine_dat$pvec >= cutoff1 & combine_dat$pvec <= cutoff2)
                mean_log_exp_bin = mean(combine_dat$target_expression[idx])
                # 
                exp_each_bin = c(exp_each_bin, mean_log_exp_bin)
            }
            
            # convert it into a dataframe and we will iterate this step for all PRTB and cell types
            exp_each_bin = as.data.frame(exp_each_bin)
            exp_each_bin$bin_order = 1:nbin
            exp_each_bin$PRTB = prtb
            exp_each_bin$CELLTYPE = CELLTYPE
            
            final_res = rbind(final_res, exp_each_bin)
        }
    }
    
    # get the ggplot2 object
    p = ggplot(data = final_res, aes(x=bin_order, y=exp_each_bin, shape=CELLTYPE, color = PRTB)) +
        geom_line( ) +
        geom_point( ) +
        theme_minimal()  +
        xlab("bins for perturbation scores") +
        ylab("expression levels of target")
    
    # is facet_wrap is set, add it to the ggplot object
    if(is.null(facet_wrap)){
        return(p)
    } else if (facet_wrap == "gene"){
        return(p + facet_wrap(~ PRTB, scales = facet_scale))
    } else if (facet_wrap == "split.by"){
        return(p + facet_wrap(~ CELLTYPE, scales = facet_scale))
    }
    
    
}



#' Single-cell heatmap for selected DE genes stratified by target expression
#' 
#' This function will generate single-cell expression heatmap for selected DE genes 
#' in cells of the same perturbation target (gRNA). This function is basically a 
#' wrapper function of the Seurat::DoHeatmap(), but with easier usage to select the 
#' cells based on given gRNA identity. Cells will be ordered in each stratification based 
#' on their perturbation scores. 
#' 
#' @export
#' 
#' @import Seurat
#' 
#' @param object a seurat object returned by RunMixscale()
#' @param assay the assay name to extract the expression level data from for plotting
#' @param slot the slot name to extract the expression level data from for plotting
#' @param labels the column name in the object's meta.data that contains the target
#' gene labels
#' @param nt.class.name the classification name of non-targeting gRNA cells
#' @param PRTB the gene name of the perturbation target to be plotted
#' @param slct_condition the name of the selected condition (e.g., cell type) to be 
#' plotted. If split.by was set to NULL during RunMixscale(), then user should leave it 
#' as the default value (i.e., "con1").
#' @param slct_features a vector of the names of the selected features (usually some DE 
#' genes) to be plotted in the heatmap. 
#' @param slct_ident if an alternative ident column instead of the "labels" column is desired,
#' users can input it here for plotting.
#' @return a ggplot2 object of the single-cell heatmap.
#' 
#' 

Mixscale_DoHeatmap = function(object = NULL, 
                               assay = "RNA", 
                               slot = "data",
                               labels = "gene", 
                               nt.class.name = "NT", 
                               PRTB = NULL, 
                               slct_condition = "con1",
                               slct_features = NULL,
                              slct_ident = NULL,
                               ...){
    # get the full object of the PRTB scores
    prtb_score <- Tool(object = object, slot = "RunMixscale")
    # check if the scores exist
    if(is.null(prtb_score)){
        stop(paste0("It seems the scores are not calculated for this object. Please check!"))
    }
    
    # get the scores for the selected PRTB and condition
    scores = prtb_score[[PRTB]][[slct_condition]][, c(1:2)]  # get the column 1 (score) and 2 (label)
    scores$cell_ID = rownames(scores)
    
    # get a subsetted seurat object for the DoHeatmap() function
    sub_obj = subset(object, cells = scores$cell_ID)

    # add the perturbation score to the meta-data:
    sub_obj = AddMetaData(sub_obj, metadata = scores)
    
    # 
    sub_obj$ident_plot2 = as.character(sub_obj[[labels]][, 1])
    
    #  get the target gene expression level
    target_expression <- GetAssayData(object = sub_obj, assay = assay, slot = slot)[PRTB, scores$cell_ID]
    target_expression = target_expression[match(colnames(sub_obj), names(target_expression))]

    idx_0 = which(target_expression == 0) 
    idx_not0 = which(target_expression != 0) 
    
    sub_obj$ident_plot2[target_expression == 0 & sub_obj[[labels]][, 1] != nt.class.name] = "expr = 0"
    sub_obj$ident_plot2[target_expression != 0 & sub_obj[[labels]][, 1] != nt.class.name] = "expr > 0"
    
    sub_obj$ident_plot2 = factor(x = as.character(sub_obj$ident_plot2), levels = c(nt.class.name, "expr > 0", "expr = 0"))
    
    weight_to_reorder <- sub_obj$pvec
    ordered.cells <- names(x = weight_to_reorder)[order(weight_to_reorder, decreasing = F)]
    
    # scale the expression data
    sub_obj <- ScaleData(object = sub_obj, features = unique(c(PRTB, slct_features)), assay = assay)
    
    ### plot all celltype together
    if(is.null(slct_ident)){
        p3 = DoHeatmap(object = sub_obj, features = unique(c(PRTB, slct_features)), label = TRUE, cells = ordered.cells, assay = assay, 
                       group.by = "ident_plot2") + ggtitle(paste0("Ordered by perturbation score"))
    } else {
        p3 = DoHeatmap(object = sub_obj, features = unique(c(PRTB, slct_features)), label = TRUE, cells = ordered.cells, assay = assay, 
                       group.by = slct_ident) + ggtitle(paste0("Ordered by perturbation score"))
    }

    return(p3)
}



#' Multi-way dotplot for multi-cell-line DE results
#' 
#' This function will generate a multi-way dotplot if the DE results produces by Run_wtDE() 
#' and get_DE_mat() contain multiple different cell-lines. The a-axis will be some selected perturbations and 
#' the y-axis will be some selected DE genes that users want to explore. Within each column, 
#' multiple dots will be displayed, with size representing the DE test DE Z-score (significance)
#' and color indicating the cell line identity. It is a good way to explore the consistency and
#' heterogeneity of DE results across perturbations/cell lines. 
#' 
#' @export
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' 
#' @param name description
#' 

DE_MultiwayPlot = function(DEG_mat = NULL, 
                           zscore_cap = 10, 
                           
                           ...){
    pruned_DEG_mat = prune_DE_mat(DEG_mat = DEG_mat, 
                                  zscore_cap = zscore_cap, 
                                  mask_target = FALSE, 
                                  p_threshold = 1, 
                                  )
    
    
}




#' Draw DE Z-score heatmap for gene signatures 
#' 
#' This function will generate a standard heatmap based on the DE Z-score heatmap.
#' Only the selected significant gene signatures will be plotted in the rows. 
#' 
#' 
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' 
#' @param obj the object produced by PCApermtest()
#' @param sig_genes the object produced by get_sig_genes() 
#' @param type the type of obj and sig_genes that are being input. The "standard" (default) indicates 
#' the results from a standard within-perturbation analysis. The "hclust" indicates the results from 
#' a hierarchical clustering analysis. The "multiCCA" indicates the results from a multiCCA analysis. 
#' @param direction to indicate whether gene signatures of both directions should be plotted, or just 
#' the up-regulated genes ("up") or down-regulated genes ("down") should be plotted. 
#' @param top_n a positive integer to indicate how many gene signatures should be plotted. If provided,
#' the top top_n genes will be selected for plotting. 
#' @param zscore_cap the cap for Z-scores in the Z-score matrix. Any Z-score that is larger than this 
#' value will be capped to this value.
#' @param labRow a boolen variable to indicate if the row names should be labelled in the heatmap or not.
#' @param output_path the directory where the heatmap will be saved. 
#' @param prefix the prefix for how the file of the heatmap should be named
#' @param height the height (in inch) for the figure
#' @param width the width (in inch) for the figure
#' 
#' @return this function returns nothing. It directly output the generated figures to the directory that 
#' a user specifies. 
#' 

DE_heatmap = function(obj = NULL, 
                      sig_genes = NULL, 
                      type = c("standard", "hclust", "multiCCA"), 
                      direction = c("both", "down", "up"), 
                      top_n = 30, 
                      zscore_cap = 15,
                      labRow = T,
                      output_path = "./",
                      prefix = "heatmap",
                      height = 15, 
                      width = 12,
                      ...){
    # adjust the parameters to make them compatitble to heatmap.2()
    if(labRow == T){
        labRow = NULL
    } else {
        labRow = NA
    }
    
    if(is.null(top_n)){
        top_n = Inf
    }
    
    # do plot accordingly
    if(type == "standard"){
        # cap the Z-scores
        tmp = as.matrix(obj$mat)
        tmp[which(tmp >= zscore_cap)] = zscore_cap
        tmp[which(tmp <= -zscore_cap)] = -zscore_cap
        
        # select the genes to plot
        if(direction == "both"){
            top_up = ifelse(top_n >= length(sig_genes$upDEGs), length(sig_genes$upDEGs), top_n)
            top_down = ifelse(top_n >= length(sig_genes$downDEGs), length(sig_genes$downDEGs), top_n)
            slct_idx = c(sig_genes$downDEGs[1:top_down], rev(sig_genes$upDEGs[1:top_up]))
            
            # color range 
            col3 = rev(brewer.pal(11, "RdBu"))
            symbreaks = T
            
        } else if (direction == "down"){
            top_down = ifelse(top_n >= length(sig_genes$downDEGs), length(sig_genes$downDEGs), top_n)
            slct_idx = sig_genes$downDEGs[1:top_down]
            
            # color range 
            col3 = rev(brewer.pal(9, "Blues"))
            symbreaks = F
            
        } else if (direction == "up"){
            top_up = ifelse(top_n >= length(sig_genes$upDEGs), length(sig_genes$upDEGs), top_n)
            slct_idx = rev(sig_genes$upDEGs[1:top_up])
            
            # color range 
            col3 = brewer.pal(9, "Reds")
            symbreaks = F
            
        }
        
        if(length(slct_idx) <= 1 | any(is.na(slct_idx))){
            return(0)
        }
        
        # the plot 
        pdf(file = paste0(output_path, "/", prefix, "_", type, ".pdf"), 
            height = height, width = width)
        heatmap.2(as.matrix(tmp[ slct_idx, , drop=F]), 
                  labRow = labRow,
                  Rowv = F, Colv  = T, 
                  dendrogram = "none", margins = c(10, 10), 
                  notecol = "black", notecex = 0.3, 
                  cexRow = 1.2, cexCol = 1.4, 
                  symbreaks = symbreaks, 
                  col = col3, trace="none", 
                  keysize = 1.5, 
                  lhei = c(1.5, 11), 
                  main = paste0(prefix))
        dev.off()
        
    } else if(type == "hclust"){
        # loop through all the clusters being identified 
        for(i in 1:length(obj$cluster_assignment)){
            # get the members of the current cluster
            cluster_idx = obj$cluster_assignment[[i]]
            
            # get the Z-score matrix 
            tmp = as.matrix( obj$mat[, cluster_idx, drop = F] )
            tmp[which(tmp >= zscore_cap)] = zscore_cap
            tmp[which(tmp <= -zscore_cap)] = -zscore_cap
            
            # select the genes to plot
            tmp_sig_genes = sig_genes[[paste0(cluster_idx, collapse = "_")]]$sig_genes
            
            if(direction == "both"){
                top_up = ifelse(top_n >= length(tmp_sig_genes$upDEGs), length(tmp_sig_genes$upDEGs), top_n)
                top_down = ifelse(top_n >= length(tmp_sig_genes$downDEGs), length(tmp_sig_genes$downDEGs), top_n)
                slct_idx = c(tmp_sig_genes$downDEGs[1:top_down], rev(tmp_sig_genes$upDEGs[1:top_up]))
                
                # color range 
                col3 = rev(brewer.pal(11, "RdBu"))
                symbreaks = T
                
            } else if (direction == "down"){
                top_down = ifelse(top_n >= length(tmp_sig_genes$downDEGs), length(tmp_sig_genes$downDEGs), top_n)
                slct_idx = tmp_sig_genes$downDEGs[1:top_down]
                
                # color range 
                col3 = rev(brewer.pal(9, "Blues"))
                symbreaks = F
                
            } else if (direction == "up"){
                top_up = ifelse(top_n >= length(tmp_sig_genes$upDEGs), length(tmp_sig_genes$upDEGs), top_n)
                slct_idx = rev(tmp_sig_genes$upDEGs[1:top_up])
                
                # color range 
                col3 = brewer.pal(9, "Reds")
                symbreaks = F
            }
            
            if(length(slct_idx) <= 1 | any(is.na(slct_idx))){
                return()
            }
            
            # the plot 
            pdf(file = paste0(output_path, "/", prefix, "_cluster_",i, "_", type, ".pdf"), 
                height = height, width = width)
            heatmap.2(as.matrix(tmp[ slct_idx, , drop=F]), 
                      labRow = labRow,
                      Rowv = F, Colv  = T, 
                      dendrogram = "none", margins = c(10, 10), 
                      notecol = "black", notecex = 0.3, 
                      cexRow = 1.2, cexCol = 1.4, 
                      symbreaks = symbreaks, 
                      col = col3, trace="none", 
                      keysize = 1.5, 
                      lhei = c(1.5, 11), 
                      main = paste0(prefix, ", ", "cluster ", i))
            dev.off()
            
        }
        
    } else if(type == "multiCCA"){
        # loop through all the programs being identified 
        for(i in 1:length(obj$program_assignment)){
            # get the members of the current cluster
            # cluster_idx = obj$program_assignment[[i]]
            
            # get the Z-score matrix 
            tmp = as.matrix( sig_genes[[paste0("Program", i)]]$perm_res$mat )
            tmp[which(tmp >= zscore_cap)] = zscore_cap
            tmp[which(tmp <= -zscore_cap)] = -zscore_cap
            
            # select the genes to plot
            tmp_sig_genes = sig_genes[[paste0("Program", i)]]$sig_genes
            
            if(direction == "both"){
                top_up = ifelse(top_n >= length(tmp_sig_genes$upDEGs), length(tmp_sig_genes$upDEGs), top_n)
                top_down = ifelse(top_n >= length(tmp_sig_genes$downDEGs), length(tmp_sig_genes$downDEGs), top_n)
                slct_idx = c(tmp_sig_genes$downDEGs[1:top_down], rev(tmp_sig_genes$upDEGs[1:top_up]))
                
                # color range 
                col3 = rev(brewer.pal(11, "RdBu"))
                symbreaks = T
                
            } else if (direction == "down"){
                top_down = ifelse(top_n >= length(tmp_sig_genes$downDEGs), length(tmp_sig_genes$downDEGs), top_n)
                slct_idx = tmp_sig_genes$downDEGs[1:top_down]
                
                # color range 
                col3 = rev(brewer.pal(9, "Blues"))
                symbreaks = F
                
            } else if (direction == "up"){
                top_up = ifelse(top_n >= length(tmp_sig_genes$upDEGs), length(tmp_sig_genes$upDEGs), top_n)
                slct_idx = rev(tmp_sig_genes$upDEGs[1:top_up])
                
                # color range 
                col3 = brewer.pal(9, "Reds")
                symbreaks = F
            }
            
            if(length(slct_idx) <= 1 | any(is.na(slct_idx))){
                return()
            }
            
            # the plot 
            pdf(file = paste0(output_path, "/", prefix, "_program_",i, "_", type, ".pdf"), 
                height = height, width = width)
            heatmap.2(as.matrix(tmp[ slct_idx, , drop=F]), 
                      labRow = labRow,
                      Rowv = F, Colv  = T, 
                      dendrogram = "none", margins = c(13, 10), 
                      notecol = "black", notecex = 0.3, 
                      cexRow = 1.2, cexCol = 1.4, 
                      symbreaks = symbreaks, 
                      col = col3, trace="none", 
                      keysize = 1.5, 
                      lhei = c(1.5, 11), 
                      main = paste0(prefix, ", ", "Program ", i))
            dev.off()
        }
    } else {
        stop("please input a valid 'type' parameter!")
    }
}


#' Dot plot for enrichment results across multiple cell types
#' 
#' This function will generate a Dot plot for the enrichment results generated by Mixscale_DEenrich().
#' It will also perform multiple testing for the P-values for all the cell types taken together 
#' (or within each cell type). 
#' 
#' @export
#' @import ggplot2
#' 
#' @param obj the list of data frames generated by Mixscale_DEenrich()
#' @param adjust.methods the method for multiple testing correction method (see 'p.adjust.methods')
#' @param adjust.split TRUR/FALSE to specify if the multiple testing correction should be done for 
#' each group separately (TRUE) or all together (FALSE)
#' @param direction a character to specify what to plot: to plot the enrichment results for 
#' the up-regulated DE genes ("up") or the down-regulated DE genes ("down").
#' @param slct_labels the selected labels (cell types) that need to be plotted
#' @param log10P_cutoff the maximum value of -log10(adjusted P-value) for plotting the size of the dot (any dot with a p-value 
#' larger than this will be set to the same size)
#' @param OR_cutoff the odds ratio (of the Fisher's exact test) cutoff for plotting the color gradient 
#' of the dot (any dot with a OR larger than this will be set to the same color)
#' 
#' @return a ggplot2 object 
#' 

DEenrich_DotPlot= function(obj, 
                           adjust.methods = "BH", 
                           # adjust.split = FALSE, 
                           direction = c("up", "down", "both"),
                           log10P_cutoff = 10, 
                           OR_cutoff = 20,
                           slct_labels = NULL, 
                           plot_title = NULL) {
    # first get the names for all the cell types that need to be plotted
    if(is.null(slct_labels)){
        slct_celltype = names(obj)
    } else{
        slct_celltype = slct_labels
    }
    
    # re-format all the results into one data frame 
    all_dat = data.frame()
    for(CELLTYPE in slct_celltype){
        # load the enrichment results
        dat = obj[[CELLTYPE]]
        dat$CELLTYPE = CELLTYPE
        
        # stack the results 
        all_dat = rbind(all_dat, dat)
        rm(dat)
    }
    
    # check if the results from single or both direction should be used
    if(direction == "up"){
        all_dat = all_dat[all_dat$direction == "upDEG", ]
    } else if (direction == "down"){
        all_dat = all_dat[all_dat$direction == "downDEG", ]
    }
    
    # Perform multiple testing correction 
    all_dat$adj_P = p.adjust(all_dat$Pval, method = adjust.methods)
    all_dat$score = -log10(all_dat$adj_P)
    all_dat$score[all_dat$score >= log10P_cutoff] = log10P_cutoff
    all_dat$OR[all_dat$OR >= OR_cutoff] = OR_cutoff
    all_dat$OR[all_dat$OR <= 1] = 1
    
    # plotting 
    all_dat$CELLTYPE = factor(all_dat$CELLTYPE, levels = slct_celltype)
    all_dat$GO_term = factor(all_dat$GO_term, levels = sort(unique(all_dat$GO_term), decreasing = T))
    
    # 
    plot_obj = ggplot(all_dat, aes(x = CELLTYPE, y = GO_term, color = score, size = OR)) + 
        geom_point(alpha = 0.6) +
        geom_point(data = subset(all_dat, score > 2), shape = "*", color = "black", size = 5, alpha = 0.6) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1), 
              strip.background = element_rect(color = "black", fill = NA), 
              plot.title = element_text(size = 15, hjust = 0.5)) +
        labs(title = plot_title,
             x = "Cell type",
             y = "Pathway gene set",
             color = expression(-log[10] ~ "(adj.P)"),
             size = "OR") +
        theme(legend.position = "right") + 
        scale_color_gradient(low="grey", high="red")
    
    return(plot_obj)
    
    
}







