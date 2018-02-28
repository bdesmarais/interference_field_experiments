###################################
#### Creating a cohort-network ####
###################################
## For New Mexico

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 02/27/2018
# Last edited on: 02/27/2018
# Last edited by: Sayali

rm(list=ls())
gc()
set.seed(312)


##
## Setup
##
library(foreign)


## Reading data with a cohort column
# In the first round, we had to manually enter cohort data

data <- read.csv("~/Dropbox/Interference_in_Field_Experiments/Analysis/coppock_replication_data/nm.replication.csv", header = TRUE)

cohort_amat <- matrix(NA, nrow(data), nrow(data))
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
    if (data$cohort[i] == data$cohort[j]){
      cohort_amat[i,j] <- 1
    } else {
      cohort_amat[i,j] <- 0
    }
  }
}
  
  
save(cohort_amat, file = "cohort_network.RData")


