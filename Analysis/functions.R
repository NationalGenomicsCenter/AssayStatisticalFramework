###### Function for bootstrapping
boot_beta_probability <- function(AP,                                           # vector of assignment probabilities
                                  n_AP,                                         # how many species to subset
                                  reps,                                         # number of bootstrap replicates
                                  threshold,                                    # what assignment threshold to consider Pr(exceed)
                                  genus_weights,                                # should we weight by genus?
                                  genera){                                      # if we are weighting by genus, provide a vector of genera (all taxa)
  
  beta_lines <- list()
  p_exceed <- vector()
  p_value <- vector()
  
  for(i in 1:reps){
    
    boot_set <- sample(AP, size = n_AP, replace = F)
    family_boot <- sample(boot_set, size = 500, replace = T)                    # sample 500 with replacement; no weighting by genus
    N_beta_boot <- fitdistr(family_boot,                                        # estimate beta distribution parameters
                            densfun = "beta", 
                            start = list(shape1 = 100, shape2 = 100))
    
    beta_lines[[i]] <- N_beta_boot                                              # save out beta parameters for plotting
    
    p_exceed[i] <- 1 - pbeta(threshold,                                         # probability of exceeding a given threshold
                             shape1 = beta_lines[[i]]$estimate[1], 
                             shape2 = beta_lines[[i]]$estimate[2])
    
    p_value[i] <- AndersonDarlingTest(AP,                                       # test if consistent with beta distribution
                                      "pbeta", 
                                      shape1 = beta_lines[[i]]$estimate[1], 
                                      shape2 = beta_lines[[i]]$estimate[2])$p.value
    
    print(paste("completed", i, "reps!"))
    
  }
  
  return(list(beta_lines, p_exceed, p_value))
}