#### SPPQ submission
## Butler and Nickerson analysis
####
# Binary Cohort network with coparty indicator
####

# Sourcing functions and data import
source("butler_functions_data.R")


#### Read the original Butler and Nicketson data
#### This is the New Mexico dataset


####
# Cohort similarity as weights with coparty indicator
####

#### Read the original Butler and Nicketson data
#### This is the New Mexico dataset

z <- data$treatment #observed treatment
y.z <- data$sb24 #observed outcome
n <- length(y.z) #number of observations
t <- length(z[z==1]) #number of treated units


load("butler_cohort_copart_network_weighted.RData")

S.ideo <- cohort_copart_amat_weighted

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
      expected.exp1.low[i] <- expected.exp1.low[i] + sum(S.ideo[i,]*zp*low_support)
      expected.exp1.high[i] <- expected.exp1.high[i] + sum(S.ideo[i,]*zp*(1-low_support))
    }
    else{
      expected.exp0[i] <- expected.exp0[i] + sum(S.ideo[i,]*zp)
      expected.exp0.low[i] <- expected.exp0.low[i] + sum(S.ideo[i,]*zp*low_support)
      expected.exp0.high[i] <- expected.exp0.high[i] + sum(S.ideo[i,]*zp*(1-low_support))
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


#### Testing and p-value calculation

# Calculate observed test statistic
exposures <- indirect.treatment(permutation = z, adj.mat = S.ideo,expected.exp0.low,expected.exp1.low,expected.exp0.high,expected.exp1.high,low_support)
test.stat <- sum((lm(y.z ~ eval(z*low_support)+eval(z*(1-low_support)) + exposures[[1]]+exposures[[2]], na.action = na.omit)$resid)^2)

pval <- numeric(nrow(parameters))

registerDoParallel(cores = ncores)

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

save(list=c("BFP.results","parameters"),file="CoppockSPPQRRresults_copartisan_cohort_weighted.RData")


