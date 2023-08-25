
b_list = list()
for(i in 1:10){
    b <- rnorm(40)
    names(b) <- LETTERS
    
    b_list[[paste0("b_", i)]] = b
}


rbo_real = sapply(X = b_list, 
                  FUN = rbo, 
                  list2 = a, 
                  p)

rbo(list1 = a, list2 = b_list[[1]], p = p)



upDEG = DEG_v1$avg_log2FC[DEG_v1$avg_log2FC < 0 ]
names(upDEG) = rownames(DEG_v1[DEG_v1$avg_log2FC < 0, ])

rbo_real = sapply(X = go_term_db2, 
                  FUN = rbo, 
                  list2 = upDEG, 
                  p = 0.995, side = "bottom", k = 300)

rbo(list1 = go_term_db2$IFNB_program1_down, 
    list2 = upDEG, 
    p = 0.995, k = 300, uneven.lengths = T, side = "bottom")  

##################
# test permutation
input_list = upDEG

set.seed(10)
n_iter = 1000


# Shuffle the input_list for n_iter times (shuffled matrix has n_iter columns)
shuffled_matrix <- replicate(n_iter, sample(input_list))

p = 0.995
k = 300
uneven.lengths = T
side = "bottom"
mid = NULL

# 
# Function to calculate rbo for each go_term against the shuffled matrix
calculate_rbo <- function(go_term) {
    apply(shuffled_matrix, 
          MARGIN = 2, 
          FUN = rbo, 
          list2 = go_term, 
          p = p, 
          k = k, 
          side = side, 
          mid = mid, 
          uneven.lengths = uneven.lengths)
}

# Apply the function to each go_term in go_term_db2
system.time(expr = {
    list_perm_vector <- lapply(go_term_db2, calculate_rbo)
})
# user  system elapsed 
# 95.301   5.219 100.528 


calculate_proportion <- function(element, list_vec) {
    mean(list_vec > element)
}

p_values <- mapply(calculate_proportion, rbo_real, list_perm_vector)

data.frame(GO_term = names(rbo_real), 
           RBO = rbo_real, 
           Pval = p_values, 
           m_GO_term = sapply(X = go_term_db2, FUN = length), 
           m_intersect = sapply(X = go_term_db, FUN = function(x, y) length(intersect(x, y)), y = names(input_list)))


###########################
#  Test on Aug 24 2023 
plist = readRDS("/Users/uqljian5/Desktop/Lab_stuffs_NYGC/Paper_information_gathering/results/gene_set_database/inhouse_database/P3_signatures_2023Jun19.rds")
plist = Reduce(c, plist)

x = plist$IFNB_program1_down
y = plist$IFNG_program1_down

list1 = 1:length(x)
names(list1) = x

list2 = 1:length(y)
names(list2) = y


rbo(list1, list2, p = 0.99, k = 300, side = "bottom")

rbo2(list1, list2, p = 0.99, k = 300, side = "bottom")



system.time(expr = {
    for(i in 1:1000){
        x = rbo(list1, list2, p = 0.99, k = 300, side = "bottom")
    }
})
# user  system elapsed 
# 3.874   0.154   4.029 



system.time(expr = {
    for(i in 1:1000){
        x = rbo2(list1, list2, p = 0.99, k = 300, side = "bottom")
    }
})
# user  system elapsed 
# 9.889   0.033   9.914 









