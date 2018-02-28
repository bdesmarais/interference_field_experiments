###################################
#### Creating a cohort-network ####
###################################
## For New Mexico

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 02/27/2018
# Last edited on: 02/28/2018
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

data <- read.csv("~/git/interference_field_experiments/Analysis/coppock_replication_data/nm.replication.csv", header = TRUE)


## Cohort network
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

save(cohort_amat, file = "~/git/interference_field_experiments/Analysis/coppock_replication_data/Extensions/cohort_network.RData")


## Cohort plus copart network
cohort_copart_amat <- matrix(NA, nrow(data), nrow(data))
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
    if (data$cohort[i] == data$cohort[j] && data$party[i] == data$party[j]){
      cohort_copart_amat[i,j] <- 1
    } else {
      cohort_copart_amat[i,j] <- 0
    }
  }
}

save(cohort_copart_amat, file = "~/git/interference_field_experiments/Analysis/coppock_replication_data/Extensions/cohort_copart_network.RData")
