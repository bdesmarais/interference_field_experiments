##########################################
#### Implementing Bowers et al method ####
##########################################

rm(list=ls())
set.seed(132)

library(doParallel)
library(foreach)
library(kSamples)
library(network)
library(permute)


############################
#### Create the network ####
############################

# 
# # At first, I will work with the numbers used in the paper
# n <- 256 #number of subjects
# t <- 128 #number of treated units
# k <- 512 #number of edges
n <- 70 #number of subjects
t <- 30 #number of treated units
k <- 50 #number of edges
p <- 0.5 #probability of receiving treatment
beta <- 1.3 #causal effect parameter
tau <- 0.7 #interference parameter
perms <- 100 #number of permutations to use in testing


# Creating vector of treatment assignment; z
z <- rep(NA, n) #a vector of observed outcomes
temp.treat <- sample(c(1:n), size=t, replace=FALSE, prob=rep(1/n, n))
for (i in 1:t){
  z[temp.treat[i]] <- 1
}
z[is.na(z)] <- 0
rm(temp.treat)


# Creating adjacency matrix; S
S <- matrix(NA, n, n)
temp.edge.index <- sample(c(1:n^2), size=k, replace=FALSE, prob=rep(1/n^2, n^2))
temp.edges <- rep(NA, n^2)
for (i in 1:k){
  temp.edges[temp.edge.index[i]] <- 1
}
temp.edges[is.na(temp.edges)] <- 0
S <- matrix(temp.edges, n, n, byrow = TRUE)
rm(temp.edge.index)
rm(temp.edges)


# Visualizing the network
network <- network(S)
plot(network)


# Creating a vector of outcomes. Normal with mean 3, variance 0.8 for now
y.z <- rnorm(n, 3, sqrt(0.8))


# scalar <- as.vector(t(z)%*%S)
# spillover <- rep(NA, n)
# for (i in 1:n){
#   spillover[i] <- beta + ((1-z[i]) * (1-beta) * exp(-tau^2 * scalar[i]))
# }
# 
# 
# Creating vector of new treatment assignment; w
# w <- rep(NA, n) #a vector of observed outcomes
# temp.treat <- sample(c(1:n), size=t, replace=FALSE, prob=rep(1/n, n))
# for (i in 1:t){
#   w[temp.treat[i]] <- 1
# }
# w[is.na(w)] <- 0
# rm(temp.treat)
# 
# 
# Creating a vector of outcomes for uniformity trial
# Normal with mean 5, variance 0.2 for now
# y.0 <- rnorm(n, 5, sqrt(0.2))
# 

############################
#### Potential outcomes ####
############################

#### Transform uniformity trial outcome into observed outcome

unif.to.z <- function(z, S, y.0, beta, tau){
  # z: observed treatment assignment
  # S: adjacency matrix
  # y.0: outcome vector for uniformity trial
  # beta: growth curve parameter
  # tau: rate of growth parameter
  
  
  scalar <- as.vector(t(z)%*%S)
  spillover <- rep(NA, n)
  
  spillover <- beta + ((1-z) * (1-beta) * exp(-tau^2 * scalar))
  
  # This is equation 4
  h.y0.z <- spillover*y.0

}


#### Transform observed outcome into uniformity trial outcome

z.to.unif <- function(z, S, y.z, beta, tau){
  # z: initial treatment assignment
  # S: adjacency matrix
  # y.z: observed outcome vector
  # beta: growth curve parameter
  # tau: rate of growth parameter

  scalar <- as.vector(t(z)%*%S)
  spillover <- rep(NA, n)
  
  # Equation (3)
  spillover <- beta + ((1-z) * (1-beta) * exp(-tau^2 * scalar))
  
  # This is equation 5
  h.yz.0 <- (1/spillover)*y.z
  return(h.yz.0)
}


#### Transform observed outcome into outcome for ANY other assignment w

z.to.w <- function(z, S, w, y.z, beta, tau){
  # z: initial treatment assignment
  # S: adjacency matrix
  # w: new treatment assignment
  # y.z: vector of outcomes for z
  # beta: growth curve parameter
  # tau: rate of growth parameter
  
  scalar.z <- as.vector(t(z)%*%S)
  scalar.w <- as.vector(t(w)%*%S)
  
  spillover.z <- rep(NA, n)
  spillover.z <- beta + ((1-z) * (1-beta) * exp(-tau^2 * scalar.z))
  
  
  spillover.w <- rep(NA, n)
  spillover.w <- beta + ((1-w) * (1-beta) * exp(-tau^2 * scalar.w))
  
  
  # Below is the actual function that transforms observed outcomes into potential outcomes
  # Equation (6)
  
  h.z.to.w <- (spillover.w / spillover.z) * y.z
  
}


#########################################
#### Testing and p-value calculation ####
#########################################


p.val <- function(z, y.z){
  
  cl <- makeCluster(4) #Setup for parallel computing
  registerDoParallel(cl)
  
  # Calculate the outcome vector after taking away the effect of treatment
  y.0 <- z.to.unif(z=z, S=S, y.z=y.z, beta=beta, tau=tau)
  
  # Calculate test statistic
  test.stat <- ks.test(y.0[z==1], y.0[z==0],
                       alternative = "two.sided")$statistic
  
  
  # Calculate a vector of test statistic using permutations
  
  results <- foreach (i = 1:perms) %dopar%{
    require(permute)
    perm.z <- z[sample(1:length(z),length(z),rep=F)] #Each time we sample a permutation of z
    perm.test.stat <- ks.test(y.0[perm.z==1], y.0[perm.z==0],
                              alternative = "two.sided")$statistic
  }
  
  stopCluster(cl)
  
  # A vector of test statistics
  all.test.stat.vals <- unlist(results)
  
  # Calculating p-value
  pval <- sum(all.test.stat.vals > test.stat)/perms
  return(pval)
}




