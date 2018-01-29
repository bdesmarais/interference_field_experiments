###########################################################
##### This code reads .RData matrix of saved p-values #####
#### Prepares a table of results for reporting purpose ####
#### Also draws the plots, in case of incomplete title #### 
###########################################################

rm(list=ls())

## This is the folder where .RData file is saved
setwd("~/Dropbox/Interference_in_Field_Experiments/Analysis/Broockman_dataverse_repository/For_cluster/Expected_exposure_analysis")

## Load the matrix of pvalues
## Loaded object will have name 'pvals'
load("pvals_broockman_demvotepct.RData")

## betas must be set to same values as relevant r code
## They can be backtraced
## All codes (Coppock, Bergan (Michigan), Broockman) search over identical grid
## (-0.5,0.5) on both betas
## dimensions of pvals will determine the vector of betas
beta1s <- seq(from=-0.5, to=0.5, length.out = dim(pvals)[1])
beta2s <- seq(from=-0.5, to=0.5, length.out = dim(pvals)[1])

## Tabulate all results
high.p.value <- max(pvals)
highest.p.indices <- which(pvals==max(pvals), arr.ind = TRUE)
direct.effect.PI <- beta1s[which(pvals==max(pvals), arr.ind = TRUE)[1]]
indirect.effect.PI <- beta2s[which(pvals==max(pvals), arr.ind = TRUE)[2]]
direct.effect.CI.high <- beta1s[max(which(pvals[,which(beta2s==indirect.effect.PI)] >= 0.05))]
direct.effect.CI.low <- beta1s[min(which(pvals[,which(beta2s==indirect.effect.PI)] >= 0.05))]
indirect.effect.CI.high <- beta2s[max(which(pvals[which(beta1s==direct.effect.PI),] >= 0.05))]
indirect.effect.CI.low <- beta2s[min(which(pvals[which(beta1s==direct.effect.PI),] >= 0.05))]
result <- rbind(high.p.value, direct.effect.PI, indirect.effect.PI,
                direct.effect.CI.low, direct.effect.CI.high,
                indirect.effect.CI.low, indirect.effect.CI.high)
result


# Save
library(fields)

pdf("pval_plot_broockman_demvotepct_blackpct_raw.pdf")
image.plot(beta1s, beta2s, pvals,
           main = "Plot of p-values: Combined network",
           xlab = "Direct effects", ylab = "Indirect effects")

# Lines for point estimate
lines(beta1s, rep(indirect.effect.PI, nrow(pvals)),
      type = "l", col = "yellow", lty = 1) #indirect

lines(rep(direct.effect.PI, nrow(pvals)), beta2s,
      type = "l", col = "yellow", lty = 1) #direct

# Lines for 95% CI
lines(beta1s, rep(indirect.effect.CI.low, nrow(pvals)),
      type = "l", col = "yellow", lty = 2) #indirect low

lines(beta1s, rep(indirect.effect.CI.high, nrow(pvals)),
      type = "l", col = "yellow", lty = 2) #indirect high

lines(rep(direct.effect.CI.high, nrow(pvals)), beta2s,
      type = "l", col = "yellow", lty = 2) #direct high

lines(rep(direct.effect.CI.low, nrow(pvals)), beta2s,
      type = "l", col = "yellow", lty = 2) #direct low

dev.off()

