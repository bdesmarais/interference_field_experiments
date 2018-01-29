# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 1 Random Assignment
# this file creates 10,000 Resticted Randomizations and outputs them, along with probabilities of assignment.

# Uncomment to set your working directory
# setwd("")

rm(list=ls())
load("Experiment_1_Adjacencies.RData")

### Replicating Randomization Protocol
# Restriction proceeds as follows:
# randomization according to the Handcoded Adjacency Matric, until we've assigned 25.
# Under the logic that if we had gotten a "bad" we would have found the bug
# This restriction kicks out any "bad" randomizations using the adj. matrix that defines adj as 1/2 mile.

set.seed(1243567)

nodes <- 1:93
sims <- 10000
Z93.block <- NULL
error.counter <- 0
i <- 0
j <- 0

while (i <sims) {
  available <- rep(1, length(nodes))
  rand <- runif(93)
  Z93 <- rep(0, length(nodes))
  while (sum(available) >0 & sum(Z93)<25){
    current <- nodes[rand==max(rand[available==1])]
    available[nodes==current] <-0
    Z93[nodes == current] <- 1
    removed <- which(A.handcoded[current, ]==1)
    available[nodes %in% removed] <-0
  }
  numberofspillovers <- A.93 %*% Z93
  checker.1.1 <- (Z93==1 &numberofspillovers >0)
  if(sum(checker.1.1)==0){
    Z93.block<- cbind(Z93.block, Z93)  
    i <- i + 1    ### Counts the Number of accepted randomizations
  }
  if(sum(checker.1.1)>0){error.counter <- error.counter + 1}
  j <- j +1    
  #### Counts the number of tries
  print(i)
}

untreatables <- c("11", "79", "80", "82", "85")  ### These are all the untreated except 70, which we don't think is parallel.
untreatables.93 <- as.numeric(untreatables) + 1  ### convert b/c fid's are numbered 0-92

Z93.block[untreatables.93, ] <- 0

Spill.block <- A.93 %*% Z93.block

condition.block <- matrix(NA, 93, sims)
condition.block[Z93.block==1 & Spill.block>0] <- 1
condition.block[Z93.block==1 & Spill.block==0] <- 2
condition.block[Z93.block==0 & Spill.block>0] <- 3
condition.block[Z93.block==0 & Spill.block==0] <- 4

probs.1.1 <- rowMeans(condition.block == 1)  ### Never Happens
probs.1.0 <- rowMeans(condition.block == 2)
probs.0.1 <- rowMeans(condition.block == 3)
probs.0.0 <- rowMeans(condition.block == 4)

which(probs.1.1>0)
###sanity check: probabilities sum to one
probs.1.1+ probs.1.0 + probs.0.1 + probs.0.0

# save(Z93.block, condition.block, probs.1.1, probs.1.0, probs.0.1, probs.0.0, 
#     file="Experiment_1_Permutation_Matrix.Rdata")


