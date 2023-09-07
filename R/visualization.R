#' @include enrichment_test.R
#' @include decomposition.R
NULL


#' This script contains multiple functions for visualization of different output in 
#' this package. The major part is for the visualization for the decomposition.R. The 
#' other part also contains the visualization for the enrichment_test.R.


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
                      height = 10, 
                      width = 6,
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
                  lhei = c(1.5, 11))
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
                      lhei = c(1.5, 11))
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
            
            # the plot 
            pdf(file = paste0(output_path, "/", prefix, "_program_",i, "_", type, ".pdf"), 
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
                      lhei = c(1.5, 11))
            dev.off()
        }
    } else {
        stop("please input a valid 'type' parameter!")
    }
}











