# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 3 Random Assignment
# this file creates 10,000 Resticted Randomizations and outputs them, along with probabilities of assignment.

# Uncomment to set your working directory
# setwd("")

rm(list=ls())

load("Experiment_3_Adjacencies.rdata")
load("Experiment_3_Past_Elections.rdata")

id <- results$precinct_id_2012

sims <- 100000
Z.block <- matrix(NA, nrow=length(id), ncol=sims)
condition.block <- matrix(NA, nrow=length(id), ncol=sims)
n_treat <- 30
i<-1

while(i<=sims){
  available <- rep(1, length(id))
  Z.avails <- rep(0, length(available))
  for (j in 1:n_treat){
    randomly_treated <- as.numeric(id %in% sample(id[available==1], 1))
    available[randomly_treated==1] <- 0
    Z.avails[randomly_treated==1] <- 1
    available[adj %*% Z.avails >0] <- 0
  }
  Z <- Z.avails
  spillovers <- adj %*% Z
  condition <- rep(NA, length(id))
  condition[Z==1] <- 3 
  condition[Z==0 & spillovers>0] <- 2
  condition[Z==0 & spillovers==0] <- 1
  Z.block[,i] <- Z
  condition.block[,i] <- condition
  i <- i+1  
 # print(i)
}


prob.3 <- rowMeans(condition.block==3, na.rm=TRUE)
prob.2 <- rowMeans(condition.block==2, na.rm=TRUE)
prob.1 <- rowMeans(condition.block==1, na.rm=TRUE)

### Choose 10,000 of these 100,000 that fits the t-ratio<1 criterion.

restricted.sims <- 10000

restricted.Z.block <- matrix(NA, nrow=length(id), ncol=restricted.sims)
restricted.condition.block <- matrix(NA, nrow=length(id), ncol=restricted.sims)

i<-1
j<-1
while(i<=restricted.sims){
available.Z <- rep(1, sims)
available.Z.index <- 1:sims
selected.Z <- which(available.Z.index %in% sample(available.Z.index[available.Z==1], 1))
available.Z[selected.Z] <-0

condition.obs <- condition.block[ , selected.Z]
weights <- rep(NA, length(id))
weights[condition.obs==3] <- 1/prob.3[condition.obs==3]
weights[condition.obs==2] <- 1/prob.2[condition.obs==2]
weights[condition.obs==1] <- 1/prob.1[condition.obs==1]

VAGub.sim <- data.frame(results, condition.obs, weights)

fit31<- lm(romney_margin ~ condition.obs==3, weights=weights, 
           data=subset(VAGub.sim, condition.obs %in% c(1,3) & prob.3<1 & prob.3>0 & prob.1<1 & prob.1>0))
t <- summary(fit31)$coefficients[2,3]

if(abs(t) < 1){
  restricted.Z.block[,i] <- Z.block[,selected.Z]
  restricted.condition.block[,i] <- condition.block[,selected.Z]
  i <- i+1
  print(paste0("Success on the ", j,'th try!'))
}
#print(t)
#print(paste0("j equals ", j))
j<-j+1
}

restricted.prob.3 <- rowMeans(restricted.condition.block==3, na.rm=TRUE)
restricted.prob.2 <- rowMeans(restricted.condition.block==2, na.rm=TRUE)
restricted.prob.1 <- rowMeans(restricted.condition.block==1, na.rm=TRUE)

# Chosen Random Assignment
set.seed(12345)
chosen.rand <- sample(1:restricted.sims, 1)
Z.obs <- restricted.Z.block[,chosen.rand]
condition.obs <- restricted.condition.block[,chosen.rand]
assigned_to_treatment <- Z.obs
weights <- rep(NA, length(id))
weights[condition.obs==3] <- 1/restricted.prob.3[condition.obs==3]
weights[condition.obs==2] <- 1/restricted.prob.2[condition.obs==2]
weights[condition.obs==1] <- 1/restricted.prob.1[condition.obs==1]

results <- within(results, {
                  Z.obs <- Z.obs
                  condition.obs <- condition.obs
                  restricted.prob.3 <- restricted.prob.3
                  restricted.prob.2 <- restricted.prob.2
                  restricted.prob.1 <- restricted.prob.1
                  assigned_to_treatment <- assigned_to_treatment
                  weights <- weights
                  })

# save(results, restricted.condition.block, restricted.Z.block, 
#     file="Experiment_3_Data_with_Random_Assignment.rdata")
