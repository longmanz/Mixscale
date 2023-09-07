#' A script to test the functions in decomposition.R 
#' 


########################
# 1. within PRTB
subset_idx = "Merged"
PATHWAY="IFNG"
DEG_idx = "allG"

######################
if(DEG_idx == "topG30"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG30_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "topG50"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG50_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "allG"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_2023Jan09.rds"), full.names = T)
}

######################
DEG_mat = readRDS(file = file_path)
celltype_list = names(DEG_mat)
PRTB_list = colnames(DEG_mat[[1]])[-1]

# 
gene_ID = DEG_mat[[1]]$gene_ID

# we will mask the PRTB target in each PRTB to 0:
for(CELLTYPE in celltype_list){
    tmp = DEG_mat[[CELLTYPE]]
    for(PRTB in PRTB_list){
        tmp[tmp$gene_ID %in% PRTB, PRTB] = 0
    }
    DEG_mat[[CELLTYPE]] = tmp
}

#####
PRTB = PRTB_list[14]   # "IRF1"

tmp=list()
for(CELLTYPE in celltype_list){
    tmp[[CELLTYPE]] = DEG_mat[[CELLTYPE]][, PRTB]
}
tmp = Reduce(cbind, tmp)
colnames(tmp) = celltype_list
rownames(tmp) = gene_ID
tmp[is.na(tmp)] = 0
tmp[tmp > 37] = 37
tmp[tmp < -37] = -37

#
rm_row_idx = which( apply(tmp, MARGIN = 1, FUN = function(x) all(abs(x) < sqrt(qchisq(0.05, df = 1, lower.tail = F)) )) )
tmp = tmp[-rm_row_idx, ]

# then check the columns if there are too many 0s
rm_idx = vector()
for(i in 1:ncol(tmp)){
    if(sum(abs(tmp[, i]) == 0 ) >= 0.8*nrow(tmp)){
        rm_idx = c(rm_idx, i)
    }
}
if(length(rm_idx) != 0){
    print(paste("Removing columns:", celltype_list[rm_idx]))
    tmp = tmp[, -rm_idx, drop = F]
}

####  Run the PCApermtest()
res = PCApermtest(mat = tmp, row_filtering_pval = 0.05, k = 1)

sig_genes = get_sig_genes(perm_obj = res, k = 1, collapse = T)

ori_db = readRDS("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/Paper_information_gathering/results/gene_set_database/inhouse_database/P1_signatures_2023Jun19.rds")

####  test visualization:
DE_heatmap(obj = res, sig_genes = sig_genes, type = "standard", direction = "both", top_n = 30)


#######################
#   2.  within celltype
p_threshold = 0.05/30000
z_threshold = sqrt(qchisq(p_threshold, df=1, lower=F))
fc_threshold = 0.2

preset_max_PC = T
dist_thres = 1 - 0.4

#
subset_idx = "Merged"
PATHWAY="IFNG"
DEG_idx = "allG"


######################
if(DEG_idx == "topG30"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG30_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "topG50"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG50_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "allG"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_2023Jan09.rds"), full.names = T)
}

# 
DEG_mat = readRDS(file = file_path)

celltype_list = names(DEG_mat)
PRTB_list = colnames(DEG_mat[[1]])[-1]

# 
gene_ID = DEG_mat[[1]]$gene_ID

# we will mask the PRTB target in each PRTB to 0:
for(CELLTYPE in celltype_list){
    tmp = DEG_mat[[CELLTYPE]]
    for(PRTB in PRTB_list){
        tmp[tmp$gene_ID %in% PRTB, PRTB] = 0
    }
    DEG_mat[[CELLTYPE]] = tmp
}

#
CELLTYPE = celltype_list[2]

tmp=DEG_mat[[CELLTYPE]]
rownames(tmp) = tmp$gene_ID
tmp=tmp[, -1]
tmp[is.na(tmp)] = 0
tmp[tmp > 37] = 37
tmp[tmp < -37] = -37

rm_idx = vector()
for(i in 1:ncol(tmp)){
    if(sum(abs(tmp[, i]) > z_threshold) < 3){
        rm_idx = c(rm_idx, i)
    }
}
if(length(rm_idx) != 0){
    tmp = tmp[, -rm_idx, drop = F]
}


######
res = DEhclust(mat = tmp)
sig_genes = get_sig_genes_DEhclust(obj = res)


#####
ori_db = readRDS("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/Paper_information_gathering/results/gene_set_database/inhouse_database/P2_signatures_2023Jun19.rds")
str(ori_db$IFNG$IFNG_BXPC3_IFNGR1_IFNGR2_JAK1_JAK2_MAFF_STAT1_up)
str(ori_db$IFNG$IFNG_BXPC3_IFNGR1_IFNGR2_JAK1_JAK2_MAFF_STAT1_down)


#####################################################################
#  3.  multiCCA testing 
get_sig = function(x, threshold = 0.05/30000){
    return(sum(x <= threshold, na.rm = T))
}

p_threshold = 0.05/30000
z_threshold = sqrt(qchisq(p_threshold, df=1, lower=F))
fc_threshold = 0.2

# 
subset_idx = "Merged"
PATHWAY="TGFB1"
DEG_idx = "allG"

######################
if(DEG_idx == "topG30"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG30_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "topG50"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_topG50_2023Jan17.rds"), full.names = T)
} else if (DEG_idx == "allG"){
    file_path = list.files(path = paste0("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/New_Ultima_data_2022Nov/Mixscape_LOOv3_2023Jan04/", subset_idx, "_DE_results/DEG_mat/"), 
                           pattern = paste0(PATHWAY, "_", subset_idx, "_DEG_mat_0.2_2023Jan09.rds"), full.names = T)
}

DEG_mat = readRDS(file = file_path)

celltype_list = names(DEG_mat)
PRTB_list = colnames(DEG_mat[[1]])[-1]

# 
gene_ID = DEG_mat[[1]]$gene_ID

# we will mask the PRTB target in each PRTB to 0:
for(CELLTYPE in celltype_list){
    tmp = DEG_mat[[CELLTYPE]]
    for(PRTB in PRTB_list){
        if(PRTB %in% tmp$gene_ID){
            tmp[tmp$gene_ID %in% PRTB, PRTB] = 0
        } 
    }
    DEG_mat[[CELLTYPE]] = tmp
}


# remove the NAs
for(i in 1:length(DEG_mat)){
    DEG_mat[[i]][is.na(DEG_mat[[i]])] = 0
}

# we also cap the z-score at +/- 30 (any value above/below that will be converted to 30 )
for(i in 1:length(DEG_mat)){
    DEG_mat[[i]][, -1][DEG_mat[[i]][, -1] >= 37] = 37
    DEG_mat[[i]][, -1][DEG_mat[[i]][, -1] <= -37] = -37
}

# remove any column with all 0s (this will cause an error when doing SVD)
DEG_mat2 = list()
for(CELLTYPE in celltype_list){
    tmp = DEG_mat[[CELLTYPE]][, -1, drop = F]
    rm_idx = vector()
    for(i in 1:ncol(tmp)){
        if(sum(abs(tmp[, i]) > sqrt(qchisq(p_threshold, df=1, lower.tail = F))) < 5){
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


# only center each column, but do not scale it.
for(i in 1:length(DEG_mat2)){
    DEG_mat2[[i]] = scale(DEG_mat2[[i]], T, F)
}

X = DEG_mat2

### test the functions
res = DEmultiCCA(X)   # works pretty well! 
sig_genes = get_sig_genes_DEmultiCCA(res)

### comparison
ori_db = readRDS("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/Paper_information_gathering/results/gene_set_database/inhouse_database/P3_signatures_2023Jun19.rds")

length(intersect(sig_genes$Program1$sig_genes$downDEGs[1:300], ori_db$TNFA$TNFA_program1_down))



###############################################
# now proceed to visualization. 










