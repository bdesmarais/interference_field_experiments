########################
#### Broockman data ####
########################
#### Analysis using individual networks of demvotepercent and blackpercent
## This code demvotepercent

rm(list=ls())
gc()
set.seed(312)

library(doParallel)
library(fields)
library(foreach)
library(magic)
library(permute)


## Data

broockman.data <- read.table("broockman_intrinsic_motivation_1.tab", sep="\t", header=TRUE)
broockman.data <- broockman.data[is.na(broockman.data$demvotepercent) == FALSE,]


## Functions

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


get.similarity <- function(x, y){
  return((2-abs(x-y))/2)
}


## Setting up parameters

z <- broockman.data$treat_out #observed treatment
y.z <- broockman.data$code_some_response_given #observed outcome
n <- length(y.z) #number of observations
t <- length(z[z==1]) #number of treated units


## Create the adjacency matrix using demvotepercent
# Code below this was used to create it

adj.mat <- NA

for (i in 1:length(unique(broockman.data$leg_state))){
  state.mat <- matrix(NA, nrow =length(broockman.data[broockman.data$leg_state==unique(broockman.data$leg_state)[i],1]), ncol = length(broockman.data[broockman.data$leg_state==unique(broockman.data$leg_state)[i],1]))
  
  for (j in 1:length(broockman.data[broockman.data$leg_state==unique(broockman.data$leg_state)[i],1])){
    
    for (k in 1:length(broockman.data[broockman.data$leg_state==unique(broockman.data$leg_state)[i],1])){
      
      state.mat[j,k] <- get.similarity(broockman.data$demvotepercent[j], broockman.data$demvotepercent[k])
      
      diag(state.mat) <- 0
    }
  }
  adj.mat <- adiag(adj.mat, state.mat)
}

adj.mat <- adj.mat[-1,-1]


##
perms <- 1000 #number of permutations to use in generating expected exposure

perm <- replicate(perms, permute.within.categories(broockman.data$leg_state,z))

expected.exp0 <- rep(0, n)
expected.exp1 <- rep(0, n)

for(p in 1:ncol(perm)){
  #zp <- permute.within.categories(data$match_category,z)
  zp <- perm[,p]
  for(i in 1:n){
    if (zp[i] == 1){
      expected.exp1[i] <- expected.exp1[i] + sum(adj.mat[i,]*zp)
    }else{
      expected.exp0[i] <- expected.exp0[i] + sum(adj.mat[i,]*zp)
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
    raw.exp[i] <- sum(adj.mat[i,]*permutation) #apply function
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


#### Testing and p-value calculation
perms.test <- 1000 #number of permutations used in testing

beta1s <- seq(from=-0.5, to=0.5, by=.025)
beta2s <- seq(from=-0.5, to=0.5, by=.025)

pvals <- matrix(NA, length(beta1s), length(beta2s))

#Use a package called Parallel
cl <- parallel::makeCluster(20) #Setup for parallel computing
registerDoParallel(cl)

start.time <- Sys.time()
pvalues.ideology <- foreach (i = 1:length(beta1s), .packages = "foreach") %dopar% {
  abc <- foreach (j = 1:length(beta2s)) %dopar% {
    
    # Calculate observed test statistic
    exposure <- indirect.treatment(permutation = z, adj.mat = adj.mat)
    test.stat <- sum((lm(y.z ~ z + exposure, na.action = na.omit)$resid)^2)
    
    # Calculate a vector of test statistic using permutations
    
    results <- foreach (k = 1:perms.test) %dopar% {
      require(permute)
      perm.z <- permute.within.categories(broockman.data$leg_state,z)
      perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = adj.mat)
      
      perm.y.0 <- y.z + (-1 * beta2s[j] * indirect.treatment(permutation = z, adj.mat = adj.mat))
      perm.y.0[z==1] <- perm.y.0[z==1] - beta1s[i]
      
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

end.time <- Sys.time()
time <- end.time - start.time
time
stopCluster(cl)

for (i in 1:length(beta1s)){
  pvals[i,] <- unlist(pvalues.ideology[i])
}


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


## Saving the p-value matrix
save(pvals, file="pvals_broockman_demvotepct.RData")
write.table(pvals, file="pvals_broockman_demvotepct.csv",
            col.names = beta2s, row.names = beta1s)


# Save
pdf("pval_plot_broockman_demvotepct.pdf")
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

dev.off()

