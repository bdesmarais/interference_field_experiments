### 80% threshold, scaled ###

rm(list=ls())
gc()
set.seed(231)


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


#### Read the original Butler and Nicketson data
#### This is the New Mexico dataset

data <- read.table("nm.replication.tab", sep="\t", header=TRUE)

z <- data$treatment #observed treatment
y.z <- data$sb24 #observed outcome
n <- length(y.z) #number of observations
t <- length(z[z==1]) #number of treated units
perms <- 10000 #number of permutations to use in generating expected exposure
perms.test <- 1000 #number of permutations used in testing


# #### Generate Similarity Scores (this code taken from CoppockJEPS_datapreparation.R)
# 
# nmhouse2008 <-read.csv("CoppockJEPS_rollcalldata.csv")
# bills <- data.frame(nmhouse2008[5:21])
# 
# ## Nominate Scores
# 
# bills_nona <- bills
# bills_nona[bills_nona==99] <- NA
# rollcalls <- rollcall(bills_nona)
# nominate_scores <- wnominate(rollcalls, polarity=c(1, 2), minvotes=10)
# dwnom_scores <- nominate_scores$legislators$coord1D

get.similarity <- function(x, y){
  return((2-abs(x-y))/2)
}

load("./Original archival/CoppockJEPS.rdata") #BD

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

# Limit to 30% most similar
threshold <- quantile(c(S.ideo),0.80)
S.ideo <- S.ideo*(S.ideo>=threshold)

#### Generate expected exposure
perm <- replicate(perms, permute.within.categories(data$match_category,z))
#perm <- Z_block

expected.exp0 <- rep(0, n)
expected.exp1 <- rep(0, n)
expected.exp0.low <- rep(0, n)
expected.exp1.low <- rep(0, n)
expected.exp0.high <- rep(0, n)
expected.exp1.high <- rep(0, n)

for(p in 1:ncol(perm)){
  #zp <- permute.within.categories(data$match_category,z)
  zp <- perm[,p]
  for(i in 1:n){
    if (zp[i] == 1){
      expected.exp1[i] <- expected.exp1[i] + sum(S.ideo[i,]*zp)
      expected.exp1.low[i] <- expected.exp1[i] + sum(S.ideo[i,]*zp*low_support)
      expected.exp1.high[i] <- expected.exp1[i] + sum(S.ideo[i,]*zp*(1-low_support))
    }
    else{
      expected.exp0[i] <- expected.exp0[i] + sum(S.ideo[i,]*zp)
      expected.exp0.low[i] <- expected.exp0[i] + sum(S.ideo[i,]*zp*low_support)
      expected.exp0.high[i] <- expected.exp0[i] + sum(S.ideo[i,]*zp*(1-low_support))
    }
  }
}

num_treat <- apply(perm,1,sum)
num_control <- apply(1-perm,1,sum)
expected.exp1 <- expected.exp1/num_treat
expected.exp0 <- expected.exp0/num_control
expected.exp1.low <- expected.exp1.low/num_treat
expected.exp0.low <- expected.exp0.low/num_control
expected.exp1.high <- expected.exp1.high/num_treat
expected.exp0.high <- expected.exp0.high/num_control


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


#### Testing and p-value calculation

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

# Calculate observed test statistic
exposures <- indirect.treatment(permutation = z, adj.mat = S.ideo,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support)
test.stat <- sum((lm(y.z ~ eval(z*low_support)+eval(z*(1-low_support)) + exposures[[1]]+exposures[[2]], na.action = na.omit)$resid)^2)

pval <- numeric(nrow(parameters))

dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages("foreach", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/" )
install.packages("doParallel", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/" )
library(foreach,lib.loc=Sys.getenv("R_LIBS_USER"))
library(doParallel,lib.loc=Sys.getenv("R_LIBS_USER"))
registerDoParallel(cores=12)

BFP.results <- foreach(i=1:nrow(parameters)) %dopar% {
  
  perm.y.0 <- z.to.unif(outcome=y.z, beta1=parameters[i,1], beta2=parameters[i,2], beta3=parameters[i,3], beta4=parameters[i,4], permutation=z, adj.mat=S.ideo,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support)
  #perm.y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = S.ideo)
  
  # Calculate a vector of test statistic using permutations
  
  beta1 <- parameters[i,1]
  beta2 <- parameters[i,2]
  beta3 <- parameters[i,3]
  beta4 <- parameters[i,4]
  
  perm.test.stats <- numeric(perms.test)
  
  for(k in 1:perms.test){
    #perm.z <- Z_block[,sample(1:10000, 1)]
    perm.z <- perm[,sample(1:perms, 1)]
    perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = S.ideo,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support)
    exposure_low <- perm.exposure[[1]]
    exposure_high <- perm.exposure[[2]]
    y.sim <- perm.y.0+beta1*perm.z*low_support + beta2*perm.z*(1-low_support) + beta3*exposure_low + beta4*exposure_high
    perm.test.stats[k] <- sum((lm(y.sim ~ eval(perm.z*low_support)+eval(perm.z*(1-low_support))+exposure_low+exposure_high , na.action = na.omit)$resid)^2)
  }
  
  # Calculating p-value
  c(sum(perm.test.stats < test.stat)/perms.test,1)
  
}

stopImplicitCluster()

save(list=c("BFP.results","parameters"),file="CoppockSPPQRRresults_threshold80_scaled.RData")
