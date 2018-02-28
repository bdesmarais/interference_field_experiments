###################################
#### Creating a cohort-network ####
###################################

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 01/29/2018
# Last edited on: 02/27/2018
# Last edited by: Sayali

rm(list=ls())
gc()
set.seed(312)


##
## Setup
##
library(foreign)


## Sourcing original data and exporting it as csv
# In the first round, we had to manually enter cohort data
# data <- read.dta("bergan.dta", convert.underscore=TRUE)
# write.csv(data, file = "bergan_data.csv")

data <- read.csv("bergan_data.csv", header = TRUE)

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


