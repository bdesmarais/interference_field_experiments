##################################################################################
#### Nickerson-Butler data; Coppock network setup & Bowers et. al. framework ####
#################################################################################
# This is the version finalized by BD
# Results in PolNet draft are based on this
# All Coppock, Bergan and Broockman analyses based on this
# To see the raw effects, or the effects in which we don't adjust for standardized exposure; skip lines 101-117

############################################
#### Replicating Coppock using our code ####
############################################

#setwd("~/Dropbox/professional/Research/Active/causalityinnetworks-agenda/Interference_in_Field_Experiments/Analysis/coppock_replication_data/") # BD
#setwd("D:/Dropbox/Interference_in_Field_Experiments/Analysis/coppock_replication_data")#SP

rm(list=ls())
gc()
set.seed(231)

library(doParallel)
library(fields)
library(foreach)
library(kSamples)
library(network)
library(permute)
library(wnominate)

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
perms.test <- 5000 #number of permutations used in testing


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


#### Generate expected and net exposure
#### This is the spillover effect model

indirect.treatment <- function(permutation, adj.mat,expected.exp0,expected.exp1){ #any treatment assignment vector and adjacency matrix can be used
 # permutation: can be the initial treatment assignment or a permutation
 #for (i in 1:n){
  # raw.exp[i] <- sum(adj.mat[i,]*permutation)
  # }
 raw.exp <- c(adj.mat%*%permutation)
 net.exp <- raw.exp - (permutation*expected.exp1 + (1-permutation)*expected.exp0)
 standard.exp <- (net.exp - mean(net.exp))/sd(net.exp) #this is the spillover or indirect effect
 return(standard.exp)
}


#### We now model the uniformity trial transformation

z.to.unif <- function(outcome, beta1, beta2, permutation, adj.mat,expected.exp0,expected.exp1){
  # outcome: vector of direct treatment outcomes
  # beta1: direct treatment effect parameter
  # beta2: indirect treatment effect parameter
  # permutation: vector of a permutation of z (can be z itself)
  # adj.mat: adjacency matrix
  
  exposure <- indirect.treatment(permutation, adj.mat)
  # This is equation 5
  h.yz.0 <- outcome - (beta1*permutation) - (beta2*exposure_high) - beta3*exposure_low
  return(h.yz.0)
}


#### Testing and p-value calculation

beta1s <- seq(from=-0.5, to=0.5, by=.025)
beta2s <- seq(from=-0.5, to=0.5, by=.025)
beta3s <- seq(from=-0.5, to=0.5, by=.025)

# beta1: direct effect of treatment
# beta2: indirect effect of high-support treated
# beta3: indirect effect of low-support treated

pvals <- matrix(NA, length(beta1s), length(beta2s))

cl <- makeCluster(2) #Setup for parallel computing
registerDoParallel(cl)

pvalues.ideology <- foreach (i = 1:length(beta1s)) %do% {
  abc <- foreach (j = 1:length(beta2s)) %do% {
 
    # Calculate the outcome vector after taking away the effect of treatment
    # y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = S.ideo)
    
    # Calculate observed test statistic
    exposure <- indirect.treatment(permutation = z, adj.mat = S.ideo)
    test.stat <- sum((lm(y.z ~ z + exposure, na.action = na.omit)$resid)^2)
    
    perm.y.0 <- y.z + (-1 * beta2s[j] * indirect.treatment(permutation = z, adj.mat = S.ideo))
    perm.y.0[z==1] <- perm.y.0[z==1] - beta1s[i]
    #perm.y.0 <- z.to.unif(outcome = y.z, beta1 = beta1s[i], beta2 = beta2s[j], permutation = z, adj.mat = S.ideo)
    
    # Calculate a vector of test statistic using permutations
    
    results <- foreach (k = 1:perms.test) %dopar% {
      require(permute)
      #perm.z <- Z_block[,sample(1:10000, 1)]
      perm.z <- perm[,sample(1:perms, 1)]
      perm.exposure <- indirect.treatment(permutation = perm.z, adj.mat = S.ideo)
      
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
  pvals[i,] <- unlist(pvalues.ideology[i])
}

pvals #rows are direct effects, columns indirect


# Results
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
save(pvals, file="pvals_coppock_replication.RData")
write.table(pvals, file="pvals_coppock_replication.csv",
            col.names = beta2s, row.names = beta1s)


# Save
pdf("pval_plot_coppock_replication.pdf")
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


