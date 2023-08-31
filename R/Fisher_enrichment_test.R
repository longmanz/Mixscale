#'
NULL

#' This script is for the conventional Fisher's exact test enrichment method (adopted from DAVID GO analysis).
#' 
#' @export
#' 
#' @param input_list the input gene list 
#' @param background the background gene list (usually the expressed genes where the 
#' input gene list is generate from, ).
#' @param go_db a list of gene-lists (GO term). It should be a list contain multiple named vector, 
#' and each vector should be a vector of multiple marker/signature genes for some biological pathway/process.
#' @param list_gene A Boolen value to indicate if the overlapping genes between the input gene list and 
#' the GO-term should be output as well.
#' @param EASE A Boolen value to indicate if the EASE correction should be applied (see 
#' https://david.ncifcrf.gov/helps/functional_annotation.html). This is useful to mitigate the 
#' small-sample inflation when the input gene list is short (e.g., < 10).
#' 
#' @return a data frame contains the enrichment test results. Each row contains the P-value and enrichment odds
#' ratio calculated from a Fisher's exact test for one GO-term in the go_db. 

go_test = function(input_list = NULL, 
                   background = NULL, 
                   go_db = NULL, 
                   list_gene = F, 
                   EASE = F){
    PT = length(background)
    
    # if go_db is a list of lists, Reduce it down to a list of vector (remove all the intermediate layers)
    while(class(go_db[[1]]) == "list"){
        go_db = Reduce(c, go_db)
    }
    
    if(list_gene == F){
        dat = matrix(nrow = length(go_db), ncol = 6)
        i = 1
        for(GO_TERM in names(go_db)){
            PH = length(intersect(go_db[[GO_TERM]], background))
            LT = length(input_list)
            
            # LH_list = intersect(input_list, go_db[[GO_TERM]])
            LH = length(intersect(input_list, go_db[[GO_TERM]]))
            
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
        dat = matrix(nrow = length(go_db), ncol = 7)
        i = 1
        for(GO_TERM in names(go_db)){
            PH = length(intersect(go_db[[GO_TERM]], background))
            LT = length(input_list)
            
            LH_list = intersect(input_list, go_db[[GO_TERM]])
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
    return(dat)
}



