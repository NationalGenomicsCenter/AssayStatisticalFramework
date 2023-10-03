################################################################################
###
### Taylor Wilcox 9/12/2023
### Interpolation of assay specificity
### Didemnidae
###
################################################################################

### Dependencies ###############################################################
require("MASS")
require("scales")
require("DescTools")

### Set working directory 
setwd("C:/Users/taylorwilcox/Box/01. taylor.wilcox Workspace/ESTCP/InterploatedSpecificity/Bootstrapping")

### Load data ##################################################################
DIME <- read.csv("DIME_summary.csv") # DIMEsummary with max AP for each species

### Beta distribution
family <- sample(DIME$AP, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(DIME$AP,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 1, shape2 = 1))

AndersonDarlingTest(DIME$AP,                                                    # test if consistent with beta distribution
                    "pbeta", 
                    shape1 = N_beta[1]$estimate[1], 
                    shape2 = N_beta[1]$estimate[2])$p.value                     # p = 0.9447

### Thresholds
1 - pbeta(0.3, shape1 = N_beta[1]$estimate[1], 
          shape2 = N_beta[1]$estimate[2])

sum(DIME$Amp > 0.3)/nrow(DIME)

1 - pbeta(0.5, shape1 = N_beta[1]$estimate[1], 
          shape2 = N_beta[1]$estimate[2])

### Bootstrapping
DIME_10_030 <- boot_beta_probability(AP = DIME$AP,                              # 100 replicates with 10 taxa, probability exceed 0.3 threshold
                                     n_AP = 10,
                                     reps = 100,
                                     threshold = 0.30)


saveRDS(DIME_10_030, file="DIME_10_030.RData")                                  # save out

DIME_50_030 <- boot_beta_probability(AP = DIME$AP,                              # 100 replicates with 50 taxa, probability exceed 0.3 threshold
                                     n_AP = 50,
                                     reps = 100,
                                     threshold = 0.30)
saveRDS(DIME_50_030, file="DIME_50_030.RData")                                  # save out

### Now looking at a range of missingness
DIME_P030 <- c()
DIME_P050 <- c()
DIME_PVALUE <- c()
sample_sizes <- c()

nrow(DIME)
for(i in 20:nrow(DIME)){                                                        # 20 to 56 taxa 
  OUT030 <- boot_beta_probability(AP = DIME$AP,                                 # 20 - 56 taxa, rep 100 times, exceed 0.3 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.30)
  
  OUT050 <- boot_beta_probability(AP = DIME$AP,                                 # 20 - 56 taxa, rep 100 times, exceed 0.5 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.50)
  
  DIME_P030 <- append(DIME_P030, OUT030[[2]])
  DIME_P050 <- append(DIME_P050, OUT050[[2]])
  DIME_PVALUE <- append(DIME_PVALUE, OUT030[[3]])
  sample_sizes <- append(sample_sizes, rep(i, 100))
}

saveRDS(DIME_P030, file="DIME_P030.RData")                                      # save out
saveRDS(DIME_P050, file="DIME_P050.RData") 
saveRDS(DIME_PVALUE, file="DIME_PVALUE.RData") 