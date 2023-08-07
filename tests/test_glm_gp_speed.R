system.time(expr = {
    fit <- glm_gp(data = count_data2[idx_for_DE, ], 
                  design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                  col_data = mat_all, 
                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                  on_disk = FALSE)
    
})
# 
# user  system elapsed 
# 744.727  23.156 769.978 

system.time(expr = {
    fit2 <- glm_gp(data = count_data2[idx_for_DE, ], 
                   design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                   col_data = mat_all, 
                   subsample = 500, 
                   size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                   on_disk = FALSE)
    
})
#
# user  system elapsed 
# 208.752  16.667 225.999 

system.time(expr = {
    fit3 <- glm_gp(data = count_data2[idx_for_DE, ], 
                   design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                   col_data = mat_all, 
                   subsample = 1000, 
                   size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                   on_disk = FALSE)
    
})
#
# user  system elapsed 
# 254.787  11.307 266.414 


#### here we remove the weight covariates from the model, and see if the overdispersion changes or not 
system.time(expr = {
    fit.nw <- glm_gp(data = count_data2[idx_for_DE, ], 
                  design = ~ 0 + cell_type + log_ct, 
                  col_data = mat_all, 
                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                  on_disk = FALSE)
    
})
#   This seems to work the best . weight is of less impact on the overdispersion estimates 
#    user  system elapsed     
## 209.242   4.441 214.282 


system.time(expr = {
    fit.nw2 <- glm_gp(data = count_data2[idx_for_DE, ], 
                     design = ~ 0 + cell_type + log_ct, 
                     col_data = mat_all, 
                     subsample = 1000, 
                     size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                     on_disk = FALSE)
})
# 
# user  system elapsed 
# 76.569  10.281  87.269 

plot(fit.nw$overdispersions, fit.nw2$overdispersions)
abline(0, 1)




#### here we remove all covariates and only keep log_ct
system.time(expr = {
    fit.min <- glm_gp(data = count_data2[idx_for_DE, ], 
                     design = ~ 1 + log_ct, 
                     col_data = mat_all, 
                     size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                     on_disk = FALSE)
    
})
#   does not work well. 
#    user  system elapsed 
## 32.891   2.816  36.564 


#####################################################
#  so now we decide to use fit.nw. what is the runtime for it?

system.time(expr = {
    fit_rough <- glm_gp(data = count_data2[idx_for_DE, ], 
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
})
# 
# user  system elapsed (with overdispersion_shrinkage = T, meaning it adds another shrinkage step)
# 351.496  20.329 373.582 
# user  system elapsed (with overdispersion_shrinkage = F, no additional step)
# 322.655   9.451 332.636 




##########################
system.time(expr = {
    fit1 <- glm_gp(data = count_data2[idx_gene_rm, ], 
                  design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                  col_data = mat_tmp, 
                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                  on_disk = FALSE)
    
})

system.time(expr = {
    fit2 <- glm_gp(data = count_data2[idx_gene_rm, ], 
                  design = ~ 0 + cell_type + weight:cell_type + log_ct, 
                  col_data = mat_tmp, 
                  overdispersion_shrinkage = F,
                  size_factors = F, # since we have log_ct in the covariates, we do not need it. 
                  on_disk = FALSE)
    
})















