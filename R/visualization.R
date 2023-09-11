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
#' 
#' @param object a seurat object returned by PRTBScoring()
#'
#' @return a ggplot2 object that contains the ridge plot.

PRTBscore_RidgePlot = function(object = NULL, 
                               labels = "gene", 
                               nt.class.name = "NT", 
                               split.by = NULL, 
                               PRTB = NULL, 
                               ...){
    # first get the full object of the PRTB scores
    prtb_score <- Tool(object = object, slot = "PRTBScoring")
    # check if the scores exist
    if(is.null(prtb_score)){
        stop(paste0("It seems the scores are not calculated for this object. Please check!"))
    }
    
    # 
    
}





#' Plot to compare the expression level of the perturbation target and the 
#' perturbation scores
#' 
#' This function will generate a connected scatterplot to compare the mean
#' expression level of the perturbation target gene within different perturbation 
#' percentile bins. After running the PRTBScoring() function, user can specify the 
#' gene name of the perturbation target and the number of bins to divide the scores 
#' into, and this function will sutomatically generate a connected scatterplot. 
#' Multiple perturbation targets and cell types are allowed. 
#' 
#' @export
#' @import ggplot2
#' @import Seurat
#' 
#' @param object a seurat object returned by PRTBScoring()
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

PRTBscore_ScatterPlot = function(object = NULL, 
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
    prtb_score <- Tool(object = object, slot = "PRTBScoring")
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
            mean_log_exp_NT = mean(combine_dat[combine_dat$gene == nt.class.name, "target_expression"])
            
            # only focus on PRTB cells and remove NT cell
            combine_dat = combine_dat[combine_dat$gene != nt.class.name, ]
            
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
        return(p + facet_wrap(~ CELLTYPE, ncol=3, scales = facet_scale))
    }
    
    
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











