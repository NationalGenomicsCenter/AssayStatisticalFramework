################################################################################
###
### Taylor Wilcox 9/12/2023
### Interpolation of assay specificity
### Cyprinidae
###
################################################################################

### Dependencies ###############################################################
require("MASS")
require("scales")
require("DescTools")

### Set working directory 
setwd("C:/Users/taylorwilcox/Box/01. taylor.wilcox Workspace/ESTCP/InterploatedSpecificity/Bootstrapping")

### Load data ##################################################################
CMCP <- read.csv("CMCP_summary.csv") # CMCP summary with max AP for each species

### Beta distribution
family <- sample(CMCP$Amp, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(CMCP$Amp,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 1, shape2 = 1))
hist(family)
AndersonDarlingTest(CMCP$Amp,                                                    # test if consistent with beta distribution
                    "pbeta", 
                    shape1 = N_beta[1]$estimate[1], 
                    shape2 = N_beta[1]$estimate[2])$p.value                     # p <0.001

### Thresholds
1 - pbeta(0.3, shape1 = N_beta[1]$estimate[1], 
               shape2 = N_beta[1]$estimate[2])

sum(CMCP$Amp > 0.3)/nrow(CMCP) # Predicted = 0.0975, observed = 0.08565602

1 - pbeta(0.5, shape1 = N_beta[1]$estimate[1], 
               shape2 = N_beta[1]$estimate[2])

sum(CMCP$Amp > 0.5)/nrow(CMCP) # Predicted = 0.00119597, observed = 0.02788722


### Bootstrapping
CMCP_10_030 <- boot_beta_probability(AP = CMCP$Amp,                              # 100 replicates with 10 taxa, probability exceed 0.3 threshold
                                     n_AP = 10,
                                     reps = 100,
                                     threshold = 0.30)


saveRDS(CMCP_10_030, file="CMCP_10_030.RData")                                  # save out

CMCP_50_030 <- boot_beta_probability(AP = CMCP$Amp,                              # 100 replicates with 50 taxa, probability exceed 0.3 threshold
                                     n_AP = 50,
                                     reps = 100,
                                     threshold = 0.30)
saveRDS(CMCP_50_030, file="CMCP_50_030.RData")                                  # save out


### Now looking at a range of missingness
CMCP_P030 <- c()
CMCP_P050 <- c()
CMCP_PVALUE <- c()
sample_sizes <- c()

nrow(CMCP)
for(i in 20:nrow(CMCP)){                                                        # 50 to 285 taxa 
  OUT030 <- boot_beta_probability(AP = CMCP$Amp,                                 # 20 - 295 taxa, rep 100 times, exceed 0.3 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.30)
  
  OUT050 <- boot_beta_probability(AP = CMCP$Amp,                                 # 20 - 285 taxa, rep 100 times, exceed 0.5 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.50)
  
  CMCP_P030 <- append(CMCP_P030, OUT030[[2]])
  CMCP_P050 <- append(CMCP_P050, OUT050[[2]])
  CMCP_PVALUE <- append(CMCP_PVALUE, OUT030[[3]])
  sample_sizes <- append(sample_sizes, rep(i, 100))
}

saveRDS(CMCP_P030, file="CMCP_P030.RData")                                      # save out
saveRDS(CMCP_P050, file="CMCP_P050.RData") 
saveRDS(CMCP_PVALUE, file="CMCP_PVALUE.RData") 
