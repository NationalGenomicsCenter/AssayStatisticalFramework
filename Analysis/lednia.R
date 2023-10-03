################################################################################
###
### Taylor Wilcox 9/12/2023
### Interpolation of assay specificity
### Nemouridae
###
################################################################################

### Dependencies ###############################################################
require("MASS")
require("scales")
require("DescTools")

### Set working directory 
setwd("C:/Users/taylorwilcox/Box/01. taylor.wilcox Workspace/ESTCP/InterploatedSpecificity/Bootstrapping")

### Load data ##################################################################
LEDN <- read.csv("LEDN_summary.csv") #LEDN_summary with max AP for each species

### Beta distribution
family <- sample(LEDN$Amp, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(LEDN$AP,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 100, shape2 = 100))

AndersonDarlingTest(LEDN$AP,                                                    # test if consistent with beta distribution
                    "pbeta", 
                    shape1 = N_beta[1]$estimate[1], 
                    shape2 = N_beta[1]$estimate[2])$p.value                     # p = 0.369

### Thresholds
1 - pbeta(0.3, shape1 = N_beta[1]$estimate[1], 
          shape2 = N_beta[1]$estimate[2])

sum(LEDN$Amp > 0.3)/nrow(LEDN) 

1 - pbeta(0.5, shape1 = N_beta[1]$estimate[1], 
          shape2 = N_beta[1]$estimate[2])

sum(LEDN$Amp > 0.5)/nrow(LEDN) 

### Bootstrapping
LEDN_10_030 <- boot_beta_probability(AP = LEDN$AP,                              # 100 replicates with 10 taxa, probability exceed 0.3 threshold
                                     n_AP = 10,
                                     reps = 100,
                                     threshold = 0.30)
saveRDS(LEDN_10_030, file="LEDN_10_030.RData")                                  # save out

LEDN_50_030 <- boot_beta_probability(AP = LEDN$AP,                              # 100 replicates with 50 taxa, probability exceed 0.3 threshold
                                     n_AP = 50,
                                     reps = 100,
                                     threshold = 0.30)
saveRDS(LEDN_50_030, file="LEDN_50_030.RData")                                  # save out

### Now looking at a range of missingness
LEDN_P030 <- c()
LEDN_P050 <- c()
LEDN_PVALUE <- c()
sample_sizes <- c()

nrow(LEDN)
for(i in 20:nrow(LEDN)){                                                        # 20 to 78 taxa 
  OUT030 <- boot_beta_probability(AP = LEDN$AP,                                 # 20 - 78 taxa, rep 100 times, exceed 0.3 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.30)
  
  OUT050 <- boot_beta_probability(AP = LEDN$AP,                                 # 20 - 78 taxa, rep 100 times, exceed 0.5 threshold
                                  n_AP = i,
                                  reps = 100,
                                  threshold = 0.50)
  
  LEDN_P030 <- append(LEDN_P030, OUT030[[2]])
  LEDN_P050 <- append(LEDN_P050, OUT050[[2]])
  LEDN_PVALUE <- append(LEDN_PVALUE, OUT030[[3]])
  sample_sizes <- append(sample_sizes, rep(i, 100))
}

saveRDS(LEDN_P030, file="LEDN_P030.RData")                                      # save out
saveRDS(LEDN_P050, file="LEDN_P050.RData") 
saveRDS(LEDN_PVALUE, file="LEDN_PVALUE.RData") 

plot(LEDN_P030 ~ sample_sizes)
plot(LEDN_P050 ~ sample_sizes)


