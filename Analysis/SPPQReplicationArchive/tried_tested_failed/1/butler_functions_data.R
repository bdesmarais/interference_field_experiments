#### Commonly used functions and data import for the Bergan analysis

# Functions
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


matrix.max <- function(x){
  # x is the matrix with respect to which you want to find the max cell
  rowmax <- which.max(apply(x,1,max))
  colmax <- which.max(x[rowmax,])
  c(rowmax,colmax)
}


get.similarity <- function(x, y){
  return((2-abs(x-y))/2)
}


#### Generate expected and net exposure
#### This is the spillover effect model

indirect.treatment <- function(permutation, adj.mat,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support){ #any treatment assignment vector and adjacency matrix can be used
  # permutation: can be the initial treatment assignment or a permutation
  #for (i in 1:n){
  # raw.exp[i] <- sum(adj.mat[i,]*permutation)
  # }
  raw.exp.low <- c(adj.mat%*%(permutation*low_support))
  raw.exp.high <- c(adj.mat%*%(permutation*(1-low_support)))
  net.exp.low <- raw.exp.low - (permutation*expected.exp1.low + (1-permutation)*expected.exp0.low)
  net.exp.high <- raw.exp.high - (permutation*expected.exp1.high + (1-permutation)*expected.exp0.high)
  standard.exp.low <- (net.exp.low - mean(net.exp.low))/sd(net.exp.low) #this is the spillover or indirect effect
  standard.exp.high <- (net.exp.high - mean(net.exp.high))/sd(net.exp.high) #this is the spillover or indirect effect
  return(list(standard.exp.low,standard.exp.high,raw.exp.low,raw.exp.high))
}


#### We now model the uniformity trial transformation

z.to.unif <- function(outcome, beta1, beta2, beta3, beta4, permutation, adj.mat,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support){
  # outcome: vector of direct treatment outcomes
  # beta1: direct treatment effect parameter
  # beta2: indirect treatment effect parameter
  # permutation: vector of a permutation of z (can be z itself)
  # adj.mat: adjacency matrix
  
  exposures <- indirect.treatment(permutation, adj.mat,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support)
  
  exposure_low <- exposures[[1]]
  exposure_high <- exposures[[2]]
  # This is equation 5
  h.yz.0 <- outcome - (beta1*permutation*low_support) -(beta2*permutation*(1-low_support)) - (beta3*exposure_low) - beta4*exposure_high
  return(h.yz.0)
}


## Importing data
# all of it
data <- read.table("nm.replication.tab", sep="\t", header=TRUE)

load("CoppockJEPS.rdata")
dwnom_scores <- CoppockJEPS$dwnom_scores
low_support <- CoppockJEPS$lowsupport

## Create an adjacency/similarity matrix using ideology
S.ideo <- matrix(NA, ncol=70, nrow=70)
for (i in 1:70){
  for (j in 1:70){
    S.ideo[i,j] <- get.similarity(dwnom_scores[i], dwnom_scores[j])
  }
}
diag(S.ideo) <- 0
S.ideo[is.na(S.ideo)==T] <- 0


## Other details
perms <- 10000 #number of permutations to use in generating expected exposure
perms.test <- 1000 #number of permutations used in testing


## Parameters

beta1s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta2s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta3s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta4s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))

parameters <- expand.grid(beta1s,beta2s,beta3s,beta4s)

# beta1: direct effect of treatment for low
# beta2: direct effect of treatment for high
# beta3: indirect effect of low-support treated
# beta4: indirect effect of high-support treated

pvals <- numeric(nrow(parameters))



