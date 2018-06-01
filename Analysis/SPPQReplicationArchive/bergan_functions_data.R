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

#### Generating expected and net exposure
#### This is the spillover effect model

indirect.treatment <- function(permutation, adj.mat, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat){
  
  #any treatment assignment vector and adjacency matrix can be used
  # permutation: can be the initial treatment assignment or a permutation
  #for (i in 1:n){
  # raw.exp[i] <- sum(adj.mat[i,]*permutation)
  # }
  
  raw.exp.dem <- c(adj.mat%*%(permutation*democrat))
  raw.exp.rep <- c(adj.mat%*%(permutation*(1-democrat)))
  net.exp.dem <- raw.exp.dem - (permutation*expected.exp1.dem + (1-permutation)*expected.exp0.dem)
  net.exp.rep <- raw.exp.rep - (permutation*expected.exp1.rep + (1-permutation)*expected.exp0.rep)
  standard.exp.dem <- (net.exp.dem - mean(net.exp.dem))/sd(net.exp.dem) #this is the spillover or indirect effect
  standard.exp.rep <- (net.exp.rep - mean(net.exp.rep))/sd(net.exp.rep) #this is the spillover or indirect effect
  return(list(standard.exp.dem, standard.exp.rep, raw.exp.dem, raw.exp.rep))
}

#### Modeling the uniformity trial transformation

z.to.unif <- function(outcome, beta1, beta2, beta3, beta4, permutation, adj.mat, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat){
  # outcome: vector of direct treatment outcomes
  # beta1: direct treatment effect parameter
  # beta2: indirect treatment effect parameter
  # permutation: vector of a permutation of z (can be z itself)
  # adj.mat: adjacency matrix
  
  exposures <- indirect.treatment(permutation, adj.mat, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat)
  
  exposure_dem <- exposures[[1]]
  exposure_rep <- exposures[[2]]
  # This is equation 5
  
  h.yz.0 <- outcome - (beta1*permutation*democrat) - (beta2*permutation*(1-democrat)) - (beta3*exposure_dem) - (beta4*exposure_rep)
  return(h.yz.0)
}


## Importing data
data <- read.dta("bergan.dta", convert.underscore=TRUE)
data <- data[1:148,]


