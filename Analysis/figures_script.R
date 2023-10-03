### Plotting

### Packages
require("MASS")
require("scales")
require("DescTools")

### Import data ###
LEDN_10_030 <- readRDS("LEDN_10_030.RData")  
LEDN_50_030 <- readRDS("LEDN_50_030.RData") 
LEDN_P030 <- readRDS("LEDN_P030.RData")                                    
LEDN_P050 <- readRDS("LEDN_P050.RData") 
LEDN_PVALUE <-readRDS("LEDN_PVALUE.RData") 

DIME_10_030 <- readRDS("DIME_10_030.RData")  
DIME_50_030 <- readRDS("DIME_50_030.RData") 
DIME_P030 <- readRDS("DIME_P030.RData")                                    
DIME_P050 <- readRDS("DIME_P050.RData") 
DIME_PVALUE <-readRDS("DIME_PVALUE.RData") 

VIR2_10_030 <- readRDS("VIR2_10_030.RData")  
VIR2_50_030 <- readRDS("VIR2_50_030.RData") 
VIR2_P030 <- readRDS("VIR2_P030.RData")                                    
VIR2_P050 <- readRDS("VIR2_P050.RData") 
VIR2_PVALUE <-readRDS("VIR2_PVALUE.RData") 

CMCP_10_030 <- readRDS("CMCP_10_030.RData")  
CMCP_50_030 <- readRDS("CMCP_50_030.RData") 
CMCP_P030 <- readRDS("CMCP_P030.RData")                                    
CMCP_P050 <- readRDS("CMCP_P050.RData") 
CMCP_PVALUE <-readRDS("CMCP_PVALUE.RData") 

LEDN <- read.csv("LEDN_summary.csv")
DIME <- read.csv("DIME_summary.csv")
VIR2 <- read.csv("VIR2_summary.csv")
CMCP <- read.csv("CMCP_summary.csv")

### FIGURE #####################################################################

par(mfrow = c(2,2))
par(mar = c(1,4,4,1))
par(las = 1)
### Lednia tumana
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 20),
     xlab = "", ylab = "Density", main = "", axes = F)
title("Lednia tumana", line = -2)
axis(2)
box()

hist(LEDN$AP, add = T, breaks = seq(0, 1, by = 0.01), col = "darkgrey", freq = F)

### Re-calculate
family <- sample(LEDN$AP, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(LEDN$AP,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 100, shape2 = 100))

lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta$estimate[1], 
            shape2 = N_beta$estimate[2]) ~ seq(0,1,0.01),
      lwd = 2, col = "darkorange")


### Didemnum perlucidum
par(mar = c(1,1,4,4))
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 20),
     xlab = "", ylab = "", main = "", axes = F)
box()
title("Didemnum perlucidum", line = -2)

hist(DIME$AP, add = T, breaks = seq(0, 1, by = 0.01), col = "darkgrey", freq = F)

### Re-calculate
family <- sample(DIME$AP, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(DIME$AP,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 100, shape2 = 100))

lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta$estimate[1], 
            shape2 = N_beta$estimate[2]) ~ seq(0,1,0.01),
      lwd = 2, col = "darkorange")

### Faxonius virilis
par(mar = c(4,4,1,1))
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 20),
     xlab = "eDNAssay AP", ylab = "Density", main = "", axes = F)
axis(1)
axis(2)
box()
title("Faxonius virilis", line = -2)

hist(VIR2$Amp, add = T, breaks = seq(0, 1, by = 0.01), col = "darkgrey", freq = F)

### Re-calculate
family <- sample(VIR2$Amp, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(VIR2$Amp,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 100, shape2 = 100))

lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta$estimate[1], 
            shape2 = N_beta$estimate[2]) ~ seq(0,1,0.01),
      lwd = 2, col = "darkorange")


### Cyprinus carpio
par(mar = c(4,1,1,4))
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 20),
     xlab = "eDNAssay AP", ylab = "", main = "", axes = F)
axis(1)
box()
title("Cyprinus carpio", line = -2)

hist(CMCP$Amp, add = T, breaks = seq(0, 1, by = 0.01), col = "darkgrey", freq = F)

### Re-calculate
family <- sample(CMCP$Amp, size = 500, replace = T)                             # sample 500 with replacement; no weighting by genus  

N_beta <- fitdistr(CMCP$Amp,                                                     # estimate beta distribution parameters
                   densfun = "beta", 
                   start = list(shape1 = 100, shape2 = 100))

lines(dbeta(seq(0,1,0.01), 
            shape1 = N_beta$estimate[1], 
            shape2 = N_beta$estimate[2]) ~ seq(0,1,0.01),
      lwd = 2, col = "darkorange")

### FIGURE #####################################################################
length(LEDN_P030)

LEDN_samplesize <- c()
for(i in 20:nrow(LEDN)){
  LEDN_samplesize <- append(LEDN_samplesize, rep(i, 100))
}

LEDN_P030_lwr <- c()
LEDN_P030_upr <- c()

for(i in 20:nrow(LEDN)){
  counter <- (i-20)*100 + 1 
  LEDN_P030_lwr[i] <- quantile(LEDN_P030[counter:(counter+99)], 0.025)
  LEDN_P030_upr[i] <- quantile(LEDN_P030[counter:(counter+99)], 0.975)
}

DIME_samplesize <- c()
for(i in 20:nrow(DIME)){
  DIME_samplesize <- append(DIME_samplesize, rep(i, 100))
}

DIME_P030_lwr <- c()
DIME_P030_upr <- c()

for(i in 20:nrow(DIME)){
  counter <- (i-20)*100 + 1 
  DIME_P030_lwr[i] <- quantile(DIME_P030[counter:(counter+99)], 0.025)
  DIME_P030_upr[i] <- quantile(DIME_P030[counter:(counter+99)], 0.975)
}

VIR2_samplesize <- c()
for(i in 20:nrow(VIR2)){
  VIR2_samplesize <- append(VIR2_samplesize, rep(i, 100))
}

VIR2_P030_lwr <- c()
VIR2_P030_upr <- c()

for(i in 20:nrow(VIR2)){
  counter <- (i-20)*100 + 1 
  VIR2_P030_lwr[i] <- quantile(VIR2_P030[counter:(counter+99)], 0.025)
  VIR2_P030_upr[i] <- quantile(VIR2_P030[counter:(counter+99)], 0.975)
}

CMCP_samplesize <- c()
for(i in 20:nrow(CMCP)){
  CMCP_samplesize <- append(CMCP_samplesize, rep(i, 100))
}

CMCP_P030_lwr <- c()
CMCP_P030_upr <- c()

for(i in 20:nrow(CMCP)){
  counter <- (i-20)*100 + 1 
  CMCP_P030_lwr[i] <- quantile(CMCP_P030[counter:(counter+99)], 0.025)
  CMCP_P030_upr[i] <- quantile(CMCP_P030[counter:(counter+99)], 0.975)
}

### Plotting

par(mfrow = c(2,2))
par(las = 1)

#LEDN
par(mar = c(4,4,1,1))
plot(LEDN_P030 ~ LEDN_samplesize,
     pch = ".")
title("Lednia tumana", line = -2)
rect(20, 0, 78, 0.0112, border = F, col =  adjustcolor( "darkred", alpha.f = 0.5))

polygon(x = c(20:nrow(LEDN),
              nrow(LEDN):20), 
        y = c(LEDN_P030_lwr[20:nrow(LEDN)],
              rev(lowess(c(LEDN_P030_upr[20:nrow(LEDN)]) ~ c(20:nrow(LEDN)))$y)),
        col =  adjustcolor( "grey", alpha.f = 0.5), border = F)

lines(lowess(c(LEDN_P030_upr[20:nrow(LEDN)]) ~ c(20:nrow(LEDN))))
lines(lowess(c(LEDN_P030_lwr[20:nrow(LEDN)]) ~ c(20:nrow(LEDN))))
lines(c(20, nrow(LEDN)),c(0.0012,0.0012), lwd = 2, lty = 3)

#DIME

plot(DIME_P030 ~ DIME_samplesize,
     pch = ".")
title("Didemnum perlucidum", line = -2)

rect(20, 0.1607, nrow(DIME), 0.1407, border = F, col =  adjustcolor( "darkred", alpha.f = 0.5))

polygon(x = c(20:nrow(DIME),
              nrow(DIME):20), 
        y = c(lowess(c(DIME_P030_lwr[20:nrow(DIME)]) ~ c(20:nrow(DIME)))$y,
              rev(lowess(c(DIME_P030_upr[20:nrow(DIME)]) ~ c(20:nrow(DIME)))$y)),
        col =  adjustcolor( "grey", alpha.f = 0.5), border = F)

lines(lowess(c(DIME_P030_upr[20:nrow(DIME)]) ~ c(20:nrow(DIME))))
lines(lowess(c(DIME_P030_lwr[20:nrow(DIME)]) ~ c(20:nrow(DIME))))
lines(c(20, nrow(DIME)),c(0.1507,0.1507), lwd = 2, lty = 3)

#VIR2
### Faxonius virilis

plot(VIR2_P030 ~ VIR2_samplesize,
     pch = ".")
title("Faxonius virilis", line = -2)

rect(20, 0.0569, nrow(VIR2), 0.0769, border = F, col =  adjustcolor( "darkred", alpha.f = 0.5))

polygon(x = c(20:nrow(VIR2),
              nrow(VIR2):20), 
        y = c(lowess(c(VIR2_P030_lwr[20:nrow(VIR2)]) ~ c(20:nrow(VIR2)))$y,
              rev(lowess(c(VIR2_P030_upr[20:nrow(VIR2)]) ~ c(20:nrow(VIR2)))$y)),
        col =  adjustcolor( "grey", alpha.f = 0.5), border = F)

lines(lowess(c(VIR2_P030_upr[20:nrow(VIR2)]) ~ c(20:nrow(VIR2))))
lines(lowess(c(VIR2_P030_lwr[20:nrow(VIR2)]) ~ c(20:nrow(VIR2))))
lines(c(20, nrow(VIR2)),c(0.0669,0.0669), lwd = 2, lty = 3)

###CMCP

plot(CMCP_P030 ~ CMCP_samplesize,
     pch = ".")
title("Cyprinus carpio", line = -2)

rect(20, 0.0975 - 0.01, nrow(CMCP), 0.0975 + 0.01, border = F, col =  adjustcolor( "darkred", alpha.f = 0.5))

polygon(x = c(20:nrow(CMCP),
              nrow(CMCP):20), 
        y = c(lowess(c(CMCP_P030_lwr[20:nrow(CMCP)]) ~ c(20:nrow(CMCP)))$y,
              rev(lowess(c(CMCP_P030_upr[20:nrow(CMCP)]) ~ c(20:nrow(CMCP)))$y)),
        col =  adjustcolor( "grey", alpha.f = 0.5), border = F)


lines(lowess(c(CMCP_P030_upr[20:nrow(CMCP)]) ~ c(20:nrow(CMCP))))
lines(lowess(c(CMCP_P030_lwr[20:nrow(CMCP)]) ~ c(20:nrow(CMCP))))
lines(c(20, nrow(CMCP)),c(0.0975,0.0957), lwd = 2, lty = 3)

###

max(c(0.0012 - LEDN_P030_lwr[20], LEDN_P030_upr[20] - 0.00112))
max(c(0.0012 - LEDN_P030_lwr[50], LEDN_P030_upr[50] - 0.00112))

max(c(0.1507 - DIME_P030_lwr[20], DIME_P030_upr[20] - 0.1507))
max(c(0.1507 - DIME_P030_lwr[50], DIME_P030_upr[50] - 0.1507))

max(c(0.0669 - VIR2_P030_lwr[20], VIR2_P030_upr[20] - 0.0669))
max(c(0.0669 - VIR2_P030_lwr[50], VIR2_P030_upr[50] - 0.0669))
