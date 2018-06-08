#### SPPQ submission
## Bergan & Cole analysis
####

# Sourcing functions and data import
source("bergan_functions_data.R")


####
# Binary network with tie if MORE THAN ONE bill co-sponsored
####

# Fixing the adjacency matrix
load("bergan_cosponsorship_network.RData")
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


#########################################
#### Testing and p-value calculation ####
#########################################

# Calculate observed test statistic

exposures <- indirect.treatment(permutation = z, adj.mat = network, expected.exp0.dem, expected.exp1.dem, expected.exp0.rep, expected.exp1.rep, democrat)

test.stat <- sum((lm(y.z ~ eval(z*democrat) + eval(z*(1-democrat)) +
                       exposures[[1]] + exposures[[2]],
                     na.action = na.omit)$resid)^2)

pval <- numeric(nrow(parameters))

registerDoParallel(cores = ncores)

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

save(list=c("BFP.results","parameters"),file="BerganSPPQRRresults_cospon_binary.RData")



