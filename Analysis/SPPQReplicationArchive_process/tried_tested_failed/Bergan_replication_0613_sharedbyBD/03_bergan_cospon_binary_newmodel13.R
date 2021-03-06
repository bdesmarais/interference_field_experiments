############################
#### Bergan (Michigan) ####
############################
# Co-sponsorship network
# Binary network with tie if MORE THAN ONE bill co-sponsored
# Reparametrized model

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 03/02/2018
# Last edited on: 03/07/2018
# Last edited by: Sayali

rm(list=ls())
gc()
set.seed(132)

#Packages
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages("iterators", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("foreign", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("foreach", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("doParallel", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
library(foreign,lib.loc=Sys.getenv("R_LIBS_USER"))
library(foreach,lib.loc=Sys.getenv("R_LIBS_USER"))
library(doParallel,lib.loc=Sys.getenv("R_LIBS_USER"))

# library(doParallel)
# library(fields)
# library(foreach)
# library(foreign)
# library(kSamples)
# library(network)
# library(permute)

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


## Importing data
data <- read.dta("bergan.dta", convert.underscore=TRUE)
data <- data[1:148,]


# Fixing the adjacency matrix
load("cosponsorship_network.RData")
network <- cosponsorship_network[rownames(cosponsorship_network)[is.na(match(rownames(cosponsorship_network),
                                                                             data$name))==FALSE],
                                 rownames(cosponsorship_network)[is.na(match(rownames(cosponsorship_network),
                                                                             data$name))==FALSE]]
rm(cosponsorship_network)
network[network == "1"] <- 0 #no tie if only one bill cosponsored
network[network != "0"] <- 1
gc()


## Cleaning it up
network <- network[-which(data$finalvote < 0), -which(data$finalvote < 0)]
data <- data[-which(data$finalvote < 0), ]


## Setting treatment and outcome vector
y.z <- data$finalvote
z <- data$anytreat
perms <- 10000 #number of permutations to use in generating expected exposure
perms.test <- 1000 #number of permutations used in testing
n <- length(z)

democrat <- data$democrat


#### Generate expected exposure
perm <- replicate(perms, permute.within.categories(data$strata,z))

expected.exp0 <- rep(0, n)
expected.exp1 <- rep(0, n)
expected.exp0.dem <- rep(0, n)
expected.exp1.dem <- rep(0, n)
expected.exp0.rep <- rep(0, n)
expected.exp1.rep <- rep(0, n)

for(p in 1:ncol(perm)){
  #zp <- permute.within.categories(data$match_category,z)
  zp <- perm[,p]
  for(i in 1:n){
    if (zp[i] == 1){
      expected.exp1[i] <- expected.exp1[i] + sum(network[i,]*zp)
      expected.exp1.dem[i] <- expected.exp1.dem[i] + sum(network[i,]*zp*democrat)
      expected.exp1.rep[i] <- expected.exp1.rep[i] + sum(network[i,]*zp*(1-democrat))
    }
    else{
      expected.exp0[i] <- expected.exp0[i] + sum(network[i,]*zp)
      expected.exp0.dem[i] <- expected.exp0.dem[i] + sum(network[i,]*zp*democrat)
      expected.exp0.rep[i] <- expected.exp0.rep[i] + sum(network[i,]*zp*(1-democrat))
    }
  }
}

num_treat <- apply(perm,1,sum)
num_control <- apply(1-perm,1,sum)
expected.exp1 <- expected.exp1/num_treat
expected.exp0 <- expected.exp0/num_control
expected.exp1.dem <- expected.exp1.dem/num_treat
expected.exp0.dem <- expected.exp0.dem/num_control
expected.exp1.rep <- expected.exp1.rep/num_treat
expected.exp0.rep <- expected.exp0.rep/num_control

for(i in 1:n){
  if(num_control[i] == 0){
    expected.exp0[i] <- 0
    expected.exp0.dem[i] <- 0
    expected.exp0.rep[i] <- 0
  } else {
    expected.exp1[i] <- 0
    expected.exp1.dem[i] <- 0
    expected.exp1.rep[i] <- 0
  }
}


#### Generate expected and net exposure
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


#### We now model the uniformity trial transformation

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



#########################################
#### Testing and p-value calculation ####
#########################################

job = 13


beta1s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta2s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta3s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))
beta4s <- c(seq(from=-0.5, to=0.5, by=.05),seq(-.025,.025,length=6))

parameters <- expand.grid(beta1s, beta2s, beta3s, beta4s)


interval <- ceiling(nrow(parameters)/16)
interval.job.lower <- interval*(job-1)+1
interval.job.upper <- min(c(interval*job,nrow(parameters)))

parameters <- parameters[interval.job.lower:interval.job.upper,]  

resultsFile <- paste("BerganSPPQRRresults_cospon_binary",job,".RData",sep="")

# beta1: direct effect of treatment for dem
# beta2: direct effect of treatment for rep
# beta3: indirect effect of dem treated
# beta4: indirect effect of rep treated

pvals <- numeric(nrow(parameters))

exposures <- indirect.treatment(permutation = z, adj.mat = network, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat)

test.stat <- sum((lm(y.z ~ eval(z*democrat) + eval(z*(1-democrat)) +
                       exposures[[1]] + exposures[[2]],
                     na.action = na.omit)$resid)^2)

pval <- numeric(nrow(parameters))

registerDoParallel(cores = 20)

BFP.results <- foreach(i=1:nrow(parameters)) %dopar% {
  
  perm.y.0 <- z.to.unif(outcome=y.z, beta1=parameters[i,1], beta2=parameters[i,2], beta3=parameters[i,3], beta4=parameters[i,4], permutation=z, adj.mat=network, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat)
  #perm.y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = network)
  
  # Calculate a vector of test statistic using permutations
  
  beta1 <- parameters[i,1]
  beta2 <- parameters[i,2]
  beta3 <- parameters[i,3]
  beta4 <- parameters[i,4]
  
  perm.test.stats <- numeric(perms.test)
  
  for(k in 1:perms.test){
    #perm.z <- Z_block[,sample(1:10000, 1)]
    perm.z <- perm[,sample(1:perms, 1)]
    perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = network, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat)
    exposure_dem <- perm.exposure[[3]]
    exposure_rep <- perm.exposure[[4]]
    y.sim <- perm.y.0 + beta1*perm.z*democrat + beta2*perm.z*(1-democrat) + beta3*exposure_dem + beta4*exposure_rep
    perm.test.stats[k] <- sum((lm(y.sim ~ eval(perm.z*democrat) + eval(perm.z*(1-democrat)) + exposure_dem + exposure_rep , na.action = na.omit)$resid)^2)
    
  }
  
  # Calculating p-value
  c(sum(perm.test.stats < test.stat)/perms.test,1)
  
}

stopImplicitCluster()

save(list=c("BFP.results","parameters"),file=resultsFile)

