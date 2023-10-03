################################################################################
###
### Taylor Wilcox 9/12/2023
### Interpolation of assay specificity
### Cambaridae
###
################################################################################

### Dependencies ###############################################################
require("MASS")
require("scales")
require("DescTools")

### Set working directory 
setwd("C:/Users/taylorwilcox/Box/01. taylor.wilcox Workspace/ESTCP/InterploatedSpecificity/Bootstrapping")

### Load data ##################################################################
VIR2 <- read.csv("VIR2_summary.csv") # VIR2 summary with max AP for each species

### Beta distribution
family <- sample(VIR2$Amp, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(VIR2$Amp,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 1, shape2 = 1))
hist(family)
AndersonDarlingTest(VIR2$Amp,                                                    # test if consistent with beta distribution
                    "pbeta", 
                    shape1 = N_beta[1]$estimate[1], 
                    shape2 = N_beta[1]$estimate[2])$p.value                     # p <0.001

### Thresholds
1 - pbeta(0.3, shape1 = N_beta[1]$estimate[1], 
               shape2 = N_beta[1]$estimate[2])

sum(VIR2$Amp > 0.3)/nrow(VIR2) # Predicted = 0.0975, observed = 0.08565602

1 - pbeta(0.5, shape1 = N_beta[1]$estimate[1], 
               shape2 = N_beta[1]$estimate[2])

sum(VIR2$Amp > 0.5)/nrow(VIR2) # Predicted = 0.00119597, observed = 0.02788722


### Bootstrapping
VIR2_10_030 <- boot_beta_probability(AP = VIR2$Amp,                              # 100 replicates with 10 taxa, probability exceed 0.3 threshold
                                     n_AP = 10,
                                     reps = 100,
                                     threshold = 0.30)


saveRDS(VIR2_10_030, file="VIR2_10_030.RData")                                  # save out

VIR2_50_030 <- boot_beta_probability(AP = VIR2$Amp,                              # 100 replicates with 50 taxa, probability exceed 0.3 threshold
                                     n_AP = 50,
                                     reps = 100,
                                     threshold = 0.30)
saveRDS(VIR2_50_030, file="VIR2_50_030.RData")                                  # save out


### Now looking at a range of missingness
VIR2_P030 <- c()
VIR2_P050 <- c()
VIR2_PVALUE <- c()
sample_sizes <- c()

nrow(VIR2)
for(i in 20:nrow(VIR2)){                                                        # 50 to 285 taxa 
  OUT030 <- boot_beta_probability(AP = VIR2$Amp,                                 # 20 - 295 taxa, rep 100 times, exceed 0.3 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.30)
  
  OUT050 <- boot_beta_probability(AP = VIR2$Amp,                                 # 20 - 285 taxa, rep 100 times, exceed 0.5 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.50)
  
  VIR2_P030 <- append(VIR2_P030, OUT030[[2]])
  VIR2_P050 <- append(VIR2_P050, OUT050[[2]])
  VIR2_PVALUE <- append(VIR2_PVALUE, OUT030[[3]])
  sample_sizes <- append(sample_sizes, rep(i, 100))
}

saveRDS(VIR2_P030, file="VIR2_P030.RData")                                      # save out
saveRDS(VIR2_P050, file="VIR2_P050.RData") 
saveRDS(VIR2_PVALUE, file="VIR2_PVALUE.RData") 
