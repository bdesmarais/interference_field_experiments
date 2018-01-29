# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 2 Random Assignment

# Uncomment to set your working directory
# setwd("")

rm(list=ls())

load("Experiment_2_Past_Elections.rdata")

id <- albany$id
musttreat <- albany$musttreat
untreatable <- albany$untreatable

sims <- 100000
Z.block <- matrix(NA, nrow=length(id), ncol=sims)
condition.block <- matrix(NA, nrow=length(id), ncol=sims)
n_treat <- 18
i<-1

while(i<=sims){
  available <- as.numeric(musttreat==0 & untreatable==0)
  Z.avails <- rep(0, length(available))
  for (j in 1:n_treat){
    randomly_treated <- as.numeric(id %in% sample(id[available==1], 1))
    available[randomly_treated==1] <- 0
    Z.avails[randomly_treated==1] <- 1
    available[adj %*% Z.avails >0] <- 0
  }
  Z <- Z.avails
  Z[musttreat==1] <-1
  spillovers <- adj %*% Z
  condition <- rep(NA, length(id))
  condition[Z==1] <- 3 
  condition[Z==0 & spillovers>0] <- 2
  condition[Z==0 & spillovers==0] <- 1
  Z.block[,i] <- Z
  condition.block[,i] <- condition
  i <- i+1  
  print(i)
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
selected.Z <- as.numeric(available.Z.index %in% sample(available.Z.index[available.Z==1], 1))
available.Z[selected.Z==1] <-0

albany.sim <- within(albany,{
condition.obs <- condition.block[,selected.Z==1]
  
weights <- rep(NA, length(id))
weights[condition.obs==3] <- 1/prob.3[condition.obs==3]
weights[condition.obs==2] <- 1/prob.2[condition.obs==2]
weights[condition.obs==1] <- 1/prob.1[condition.obs==1]
})

fit31<- lm(percentageuniverse ~ condition.obs==3, weights=weights, 
           data=subset(albany.sim, condition.obs %in% c(1,3) & prob.3<1 & prob.3>0 & prob.1<1 & prob.1>0))
t <- summary(fit31)$coefficients[2,3]
if(abs(t) < 1){
  restricted.Z.block[,i] <- Z.block[,selected.Z==1]
  restricted.condition.block[,i] <- condition.block[,selected.Z==1]
  i <- i+1
  print(paste0("Success on the ", j,'th try!'))
}
print(t)
print(paste0("j equals ", j))
j<-j+1
}


#save(restricted.condition.block, restricted.Z.block, condition.block, Z.block, 
#     file="Experiment_2_Permutation_Matrix.rdata")

load("Experiment_2_Permutation_Matrix.rdata")

restricted.prob.3 <- rowMeans(restricted.condition.block==3, na.rm=TRUE)
restricted.prob.2 <- rowMeans(restricted.condition.block==2, na.rm=TRUE)
restricted.prob.1 <- rowMeans(restricted.condition.block==1, na.rm=TRUE)

# Actual Assigment
set.seed(12345)
chosen.rand <- sample(1:restricted.sims, 1)
Z.obs <- restricted.Z.block[,chosen.rand]
assigned_to_treatment <- Z.obs
treatment_districts <- subset(data.frame(albany, assigned_to_treatment), assigned_to_treatment==1)
# write.csv(treatment_districts, "albanytreatmentdistricts.csv")


