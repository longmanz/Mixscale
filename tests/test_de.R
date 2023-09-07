library(Seurat)
library(SeuratObject)
library(glmGamPoi)
source("/Users/uqljian5/Documents/github_repo/Perturbation_Scoring/R/get_fold_change.R")
source("/Users/uqljian5/Documents/github_repo/Perturbation_Scoring/R/perturbation_scoring.R")
source("/Users/uqljian5/Documents/github_repo/Perturbation_Scoring/R/scoring_de.R")
source("/Users/uqljian5/Documents/github_repo/Perturbation_Scoring/R/glm_gp_disp_only.R")


test = readRDS("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/Paper_information_gathering/results/gRNA_efficiency/plot_2023Jun05/IFNG_Batch2_LOOv3_Parse_2023Jan04.rds")
test = test[[1]]
test = subset(test, subset = gene %in% c("NT", "IFNGR1", "ZC3H3", "TSC22D1"))

old_score = Tool(test, slot = "RunMixscape_LOOv3")
saveRDS(test, file = "/Users/uqljian5/Documents/github_repo/test_data/subset_IFNG_Batch2_Parse_2023Jan04.rds")
saveRDS(old_score, file = "/Users/uqljian5/Documents/github_repo/test_data/old_score_IFNG_Batch2_Parse_2023Jan04.rds")


object = readRDS(file = "/Users/uqljian5/Documents/github_repo/test_data/subset_IFNG_Batch2_Parse_2023Jan04.rds")
DefaultAssay(object) = "RNA"
object[['PRTB']] = NULL

object <- CalcPerturbSig(
    object = object, 
    assay = "RNA", 
    slot = "data", 
    gd.class ="gene", 
    nt.cell.class = "NT", 
    reduction = "pca", 
    ndims = 40, 
    num.neighbors = 20, 
    new.assay.name = "PRTB", 
    split.by = "cell_type")


object = PRTBScoring(
    object = object, 
    assay = "PRTB", 
    slot = "scale.data", 
    labels = "gene", 
    nt.class.name = "NT", 
    min.de.genes = 5, 
    split.by = "cell_type", 
    logfc.threshold = 0.2,
    de.assay = "RNA",
    max.de.genes = 100, 
    prtb.type = "P", new.class.name = "mixscape_v1", fine.mode = F, fine.mode.labels = "NT", harmonize = T, seed = 1)


new_score = Tool(object, slot = "PRTBScoring")
new_score2 = Tool(test2, slot = "PRTBScoring")



plot(x = old_score$ZC3H3$A549$`EMC3-AS1`, 
     y = new_score$ZC3H3$A549$`EMC3-AS1`)
abline(0, 1)


plot(x = new_score$IFNGR1$BXPC3$pvec, 
     y = new_score2$IFNGR1$BXPC3$pvec)
abline(0, 1)


#################
x = FoldChange_new(object = GetAssayData(object = test1, slot = "data"),
               cells.1 = colnames(test1)[test1$gene == "NT"],
               cells.2 = colnames(test1)[test1$gene == "IFNGR1"],
               mean.fxn = function(x) log(x = rowMeans(x = expm1(x = x)) + 0.1, base = 2),
               fc.name = "avg_log2FC",
               features = rownames(x = test1) ) 


res = scoringDE(object = object, assay = "PRTB", slot = "data", labels = "gene", 
                nt.class.name = "NT", verbose = FALSE, 
                split.by = "cell_type",
                total_ct_labels = "nCount_RNA",  
                logfc.threshold = 0.2, 
                pseudocount.use = 0.1, 
                base = 2,
                min.pct = 0.1, 
                min.cells = 10)


