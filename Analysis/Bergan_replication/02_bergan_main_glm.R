######################################################################
#### Analyze the Bergan (New Hampshire) data for spillover effect ####
######################################################################
## With glm residualsuals used for calculating test statistic

rm(list=ls())
set.seed(132)

setwd("D:/Dropbox/Interference_in_Field_Experiments/Analysis/Bergan_replication") #SP
library(doParallel)
library(fields)
library(foreach)
library(foreign)
library(kSamples)
library(network)
library(permute)

par(mfrow=c(1,1))


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
data <- data[data$finalvote >= 0,]
data <- data[order(data$newtreat),]
data <- data[1:143,]


## Setting treatment, vector and wnominate vector (for network)
y.z <- data$finalvote
z <- data$anytreat
dwnom_scores <- data$wnominate
perms <- 10000 #number of permutations to use in generating expected exposure
perms.test <- 1000 #number of permutations used in testing
n <- length(z)


get.similarity <- function(x, y){
  return((2-abs(x-y))/2)
}


## Create an adjacency/similarity matrix using ideology
S.ideo <- matrix(NA, ncol=nrow(data), nrow=nrow(data))
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
    S.ideo[i,j] <- get.similarity(dwnom_scores[i], dwnom_scores[j])
  }
}
diag(S.ideo) <- 0
S.ideo[is.na(S.ideo)==T] <- 0


#### Generate expected exposure
perm <- replicate(perms, z[sample(1:length(z),length(z),rep=F)])

expected.exp0 <- rep(0, n)
expected.exp1 <- rep(0, n)

for(p in 1:ncol(perm)){
  zp <- perm[,p]
  for(i in 1:n){
    if (zp[i] == 1){
      expected.exp1[i] <- expected.exp1[i] + sum(S.ideo[i,-i]*zp[-i])
    }
    else{
      expected.exp0[i] <- expected.exp0[i] + sum(S.ideo[i,]*zp)
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
  standard.exp <- (net.exp - mean(net.exp))/sd(net.exp) #this is the spillover or indirect effect
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

beta1s <- seq(from=-1.5, to=1.5, by=.2)
beta2s <- seq(from=-1.5, to=1.5, by=.2)

pvals <- matrix(NA, length(beta1s), length(beta2s))

cl <- makeCluster(4) #Setup for parallel computing
registerDoParallel(cl)

pvalues.ideology <- foreach (i = 1:length(beta1s)) %do% {
  abc <- foreach (j = 1:length(beta2s)) %do% {
    
    # Calculate the outcome vector after taking away the effect of treatment
    y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = S.ideo)
    
    # Calculate observed test statistic
    exposure <- indirect.treatment(permutation = z, adj.mat = S.ideo)
    test.stat <- sum((glm(y.z ~ z + exposure, na.action = na.omit,
                          family = "binomial")$residuals)^2)
    
    # Calculate a vector of test statistic using permutations
    
    results <- foreach (k = 1:perms.test) %dopar% {
      require(permute)
      perm.z <- permute.within.categories(data$strata,z)
      perm.y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = S.ideo)
      perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = S.ideo)
      
      y.sim <- plogis(perm.y.0 + beta1s[i]*perm.z + beta2s[j]*perm.exposure)
      perm.test.stat <- sum((glm(y.sim ~ perm.z + perm.exposure, na.action = na.omit,
                                family = "binomial")$residuals)^2)
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
  pvals[i,] <- unlist(pvalues.ideology[i])
}

pvals


## Creating a plot
image.plot(beta1s, beta2s, pvals,
           main = "Plot of p-values",
           xlab = "Direct effects", ylab = "Indirect effects")


pdf("pvalues_figure_main_glm.pdf")
image.plot(beta1s, beta2s, pvals,
          main = "Plot of p-values: glm",
          xlab = "Direct effects", ylab = "Indirect effects")
dev.off()

