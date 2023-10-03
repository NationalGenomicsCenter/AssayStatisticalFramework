################################################################################
###
### Taylor Wilcox 9/22/2023
### Interpolation of assay specificity
### Cyprinidae - alternative weightings
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
Genera_list <- read.csv("CMCP_genera.csv") # list of all species with a column for genus

### Format data for weighting
species_list_genus_counts <- aggregate(Genera_list$species, by = list(Genera_list$genus), FUN = length)
AP_list_genus_counts <- aggregate(CMCP$Taxon, by = list(CMCP$Genus), FUN = length)
species_list_genus_counts$weight <- species_list_genus_counts$x/sum(species_list_genus_counts$x)

combined_lists <- merge(species_list_genus_counts, CMCP, by.x = "Group.1", by.y = "Genus")

### Beta distribution
family1 <- sample(combined_lists$Amp, size = 500, replace = T)                                   # sample 500 with replacement; weighted by known species  

N_beta1 <- fitdistr(family1,                                                     # estimate beta distribution parameters
                    densfun = "beta", 
                    start = list(shape1 = 1, shape2 = 1))

family2 <- sample(combined_lists$Amp, size = 500, replace = T, 
                 prob = combined_lists$weight)                                   # sample 500 with replacement; weighted by known species  

N_beta2 <- fitdistr(family2,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 1, shape2 = 1))

### Figure
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 15),
     xlab = "eDNAssay AP", ylab = "", main = "", axes = F)
axis(1)
box()
title("Cyprinus carpio", line = -2)

hist(CMCP$Amp, add = T, breaks = seq(0, 1, by = 0.01), col = "darkgrey", freq = F)

### Re-calculate
lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta1$estimate[1], 
            shape2 = N_beta1$estimate[2]) ~ seq(0,1,0.01),
      lwd = 5, col = "darkorange")

lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta2$estimate[1], 
            shape2 = N_beta2$estimate[2]) ~ seq(0,1,0.01),
      lwd = 3, col = "blue", lty = 3)

legend("right", lwd = 3, lty = c(1,3), col = c("darkorange","blue"), 
       legend = c("Unweighted", "weighted"))
