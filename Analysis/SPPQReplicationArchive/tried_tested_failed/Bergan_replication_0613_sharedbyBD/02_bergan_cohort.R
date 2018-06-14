##################################################################
#### Analyze the Bergan (Michigan) data for spillover effect ####
##################################################################
# Cohort network
# Skip lines 70-84 for raw exposure results
# Base model replicated from Coppock

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 01/28/2018
# Last edited on: 03/02/2018
# Last edited by: Sayali

rm(list=ls())
gc()
set.seed(132)


library(doParallel)
library(fields)
library(foreach)
library(foreign)
library(kSamples)
library(network)
library(permute)


permute.within.categories <- function(categories,z){
  ucategories <- unique(categories)
  perm.z <- rep(NA,length(z))
  for(c in ucategories){
    z.c <- z[which(categories==c)]
    perm.z.c <- sample(z.c,length(z.c),rep=F)
    perm.z[which(categories==c)] <- perm.z.c
  }
  perm.z
}


## Importing data
data <- read.dta("bergan.dta", convert.underscore=TRUE)
data <- data[1:148,]

# Fixing the adjacency matrix
load("w_cohort_network.RData")
network <- w_cohort_amat
rm(w_cohort_amat)
gc()

## Cleaning it up
network <- network[-which(data$finalvote < 0), -which(data$finalvote < 0)]
data <- data[-which(data$finalvote < 0), ]


## Setting treatment and outcome vector
y.z <- data$finalvote
z <- data$anytreat
perms <- 1000 #number of permutations to use in generating expected exposure
perms.test <- 500 #number of permutations used in testing
n <- length(z)

beta1s <- seq(from=-0.5, to=0.5, by=.025)
beta2s <- seq(from=-0.5, to=0.5, by=.025)


#### Generate expected exposure
perm <- replicate(perms, permute.within.categories(data$strata,z))

expected.exp0 <- rep(0, n)
expected.exp1 <- rep(0, n)

for(p in 1:ncol(perm)){
  zp <- perm[,p]
  for(i in 1:n){
    if (zp[i] == 1){
      expected.exp1[i] <- expected.exp1[i] + sum(network[i,]*zp)
    }
    else{
      expected.exp0[i] <- expected.exp0[i] + sum(network[i,]*zp)
    }
  }
}
num_treat <- apply(perm,1,sum)
num_control <- apply(1-perm,1,sum)
expected.exp1 <- expected.exp1/num_treat
expected.exp0 <- expected.exp0/num_control


#### Generate expected and net exposure
#### This is the spillover effect model

indirect.treatment <- function(permutation, adj.mat){ #any treatment assignment vector and adjacency matrix can be used
  # permutation: can be the initial treatment assignment or a permutation
  raw.exp <- rep(NA, n)
  for (i in 1:n){
    raw.exp[i] <- sum(adj.mat[i,]*permutation)
  }
  
  net.exp <- raw.exp - (permutation*expected.exp1 + (1-permutation)*expected.exp0)
  standard.exp <- (net.exp - mean(net.exp, na.rm = TRUE))/sd(net.exp, na.rm = TRUE) #this is the spillover or indirect effect
  return(standard.exp)
}


#### We now model the uniformity trial transformation

z.to.unif <- function(outcome, beta1, beta2, permutation, adj.mat){
  # outcome: vector of direct treatment outcomes
  # beta1: direct treatment effect parameter
  # beta2: indirect treatment effect parameter
  # permutation: vector of a permutation of z (can be z itself)
  # adj.mat: adjacency matrix
  
  exposure <- indirect.treatment(permutation, adj.mat)
  # This is equation 5
  h.yz.0 <- outcome - (beta1*permutation) - (beta2*exposure)
  return(h.yz.0)
}


#########################################
#### Testing and p-value calculation ####
#########################################

pvals <- matrix(NA, length(beta1s), length(beta2s))

cl <- makeCluster(4) #Setup for parallel computing
registerDoParallel(cl)

pvaluenetwork <- foreach (i = 1:length(beta1s)) %do% {
  abc <- foreach (j = 1:length(beta2s)) %do% {
    
    # Calculate observed test statistic
    exposure <- indirect.treatment(permutation = z, adj.mat = network)
    test.stat <- sum((lm(y.z ~ z + exposure, na.action = na.omit)$resid)^2)
    
    perm.y.0 <- y.z + (-1 * beta2s[j] * indirect.treatment(permutation = z, adj.mat = network))
    perm.y.0[z==1] <- perm.y.0[z==1] - beta1s[i]
    
    # Calculate a vector of test statistic using permutations
    
    results <- foreach (k = 1:perms.test) %dopar% {
      require(permute)
      perm.z <- permute.within.categories(data$strata,z)
      perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = network)
      
      y.sim <- perm.y.0 + beta1s[i]*perm.z + beta2s[j]*perm.exposure
      perm.test.stat <- sum((lm(y.sim ~ perm.z + perm.exposure, na.action = na.omit)$resid)^2)
    }
    
    # A vector of test statistics
    all.test.stat.vals <- as.numeric(unlist(results))
    
    # Calculating p-value
    pval <- sum(all.test.stat.vals < test.stat)/perms.test
  }
  as.numeric(unlist(abc))
}

stopCluster(cl)

for (i in 1:length(beta1s)){
  pvals[i,] <- unlist(pvaluenetwork[i])
}

#pvals #rows are direct effects, columns indirect


# Saving results
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


## Creating a plot
image.plot(beta1s, beta2s, pvals,
           main = "Plot of p-values",
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


# ## Saving the p-value matrix
# save(pvals, file="pvals_bergan_cohort_raw.RData")
# write.table(pvals, file="pvals_bergan_cohort_raw.csv",
#             col.names = beta2s, row.names = beta1s)
# 
# 
# # Save
# pdf("pval_plot_bergan_cohort_raw.pdf")
# image.plot(beta1s, beta2s, pvals,
#            main = "Plot of p-values",
#            xlab = "Direct effects", ylab = "Indirect effects")
# 
# # Lines for point estimate
# lines(beta1s, rep(indirect.effect.PI, nrow(pvals)),
#       type = "l", col = "yellow", lty = 1) #indirect
# 
# lines(rep(direct.effect.PI, nrow(pvals)), beta2s,
#       type = "l", col = "yellow", lty = 1) #direct
# 
# # Lines for 95% CI
# lines(beta1s, rep(indirect.effect.CI.low, nrow(pvals)),
#       type = "l", col = "yellow", lty = 2) #indirect low
# 
# lines(beta1s, rep(indirect.effect.CI.high, nrow(pvals)),
#       type = "l", col = "yellow", lty = 2) #indirect high
# 
# lines(rep(direct.effect.CI.high, nrow(pvals)), beta2s,
#       type = "l", col = "yellow", lty = 2) #direct high
# 
# lines(rep(direct.effect.CI.low, nrow(pvals)), beta2s,
#       type = "l", col = "yellow", lty = 2) #direct low
# 
# dev.off()



