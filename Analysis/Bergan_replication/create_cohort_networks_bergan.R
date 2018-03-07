###################################
#### Creating a cohort-network ####
###################################

# Authors: Sayali Phadke, Bruce Desmarais
# Created on: 01/29/2018
# Last edited on: 03/02/2018
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
data <- data[1:148, ]


####
#### Creating and saving networks
####

## Simple cohort network
# cohort_amat <- matrix(NA, nrow(data), nrow(data))
# for (i in 1:nrow(data)){
#   for (j in 1:nrow(data)){
#     if (data$cohort[i] == data$cohort[j]){
#       cohort_amat[i,j] <- 1
#     } else {
#       cohort_amat[i,j] <- 0
#     }
#   }
# }
#   
# save(cohort_amat, file = "cohort_network.RData")




## Cohort and coparty network
cohort_copart_amat <- matrix(NA, nrow(data), nrow(data))
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
    if (data$cohort[i] == data$cohort[j] && data$democrat[i] == data$democrat[j] && data$senate[i] == data$senate[j]){
      cohort_copart_amat[i,j] <- 1
    } else {
      cohort_copart_amat[i,j] <- 0
    }
  }
}

save(cohort_copart_amat, file = "cohort_copart_network.RData")




## Weighted cohort network
# w_cohort_amat <- matrix(NA, nrow(data), nrow(data))
# for (i in 1:nrow(data)){
#   for (j in 1:nrow(data)){
#       w_cohort_amat[i,j] <- 1/(1+(abs(data$cohort[i] - data$cohort[j])))
#   }
# }
# 
# save(w_cohort_amat, file = "w_cohort_network.RData")




## Weighted cohort and coparty network
w_cohort_copart_amat <- matrix(NA, nrow(data), nrow(data))
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
    if (data$democrat[i] == data$democrat[j] && data$senate[i] == data$senate[j]){
    w_cohort_copart_amat[i,j] <- 1/(1+(abs(data$cohort[i] - data$cohort[j])))
    } else {
      w_cohort_copart_amat[i,j] <- 0
    }
  }
}

save(w_cohort_copart_amat, file = "w_cohort_copart_network.RData")
