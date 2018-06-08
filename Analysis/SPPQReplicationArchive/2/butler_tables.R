####
#### Creating tables for Butler and Nickerson analysis using results
####

rm(list=ls())
gc()


## Create a function to evaluate the point estimate and confidence intervals using the grid of p-values
BFP.results.summary <- function(parameters, p.values, level){
  threshold <- 1-level
  estimate <- as.numeric(parameters[which.max(p.values),])
  CIs <- NULL
  for(p in 1:ncol(parameters)){
    select.col <- 0
    for(c in (1:ncol(parameters))[-p]){
      select.col <- select.col + 1*(parameters[,c]==estimate[c])
    }
    select.col <- which(select.col==(ncol(parameters)-1))
    parameters.p <- parameters[select.col,p]
    p.values.p <- p.values[select.col]
    parameters.p <- parameters.p[which(p.values.p>threshold)]
    CIs <- rbind(CIs,c(min(parameters.p),max(parameters.p)))
  }
  list(estimate,CIs)
}



## Read in results of all four Butler-Nickerson analyses
# Binary committee with copartisonship
load("CoppockSPPQRRresults_copartisan_committee_binary.RData")
copartisan_committee_binary_results <- do.call('rbind',BFP.results)[,1]


# Weighted committee with copartisonship
load("CoppockSPPQRRresults_copartisan_committee_weighted.RData")
copartisan_committee_weighted_results <- do.call('rbind',BFP.results)[,1]


# Binary cohort with copartisanship
load("CoppockSPPQRRresults_copartisan_cohort_binary.RData")
copartisan_cohort_binary_results <- do.call('rbind',BFP.results)[,1]


# Weighted cohort with copartisanship
load("CoppockSPPQRRresults_copartisan_cohort_weighted.RData")
copartisan_cohort_weighted_results <- do.call('rbind',BFP.results)[,1]


## Generate the point estimates and 90 & 95% CIs for each analysis
# Binary cohort with copartisanship
summary.cohort_binary.95 <- BFP.results.summary(parameters, copartisan_cohort_binary_results, level = 0.95)
summary.cohort_binary.9 <- BFP.results.summary(parameters, copartisan_cohort_binary_results, level = 0.9)

cohort.binary.table <- cbind(summary.cohort_binary.95[[1]], summary.cohort_binary.95[[2]], summary.cohort_binary.9[[2]])
cohort.binary <- xtable(cohort.binary.table)
print(cohort.binary, file = "butler_cohort_binary_table4_bottomleft.txt")


# Weighted cohort with copartisanship
summary.cohort_weighted.95 <- BFP.results.summary(parameters, copartisan_cohort_weighted_results, level = 0.95)
summary.cohort_weighted.9 <- BFP.results.summary(parameters, copartisan_cohort_weighted_results,level = 0.9)

cohort.weighted.table <- cbind(summary.cohort_weighted.95[[1]], summary.cohort_weighted.95[[2]], summary.cohort_weighted.9[[2]])
cohort.weighted<- xtable(cohort.weighted.table)
print(cohort.weighted, file = "butler_cohort_weighted_table4_bottomright.txt")


# Binary committee with copartisanship
summary.committee_binary.95 <- BFP.results.summary(parameters, copartisan_committee_binary_results, level = 0.95)
summary.committee_binary.9 <- BFP.results.summary(parameters, copartisan_committee_binary_results,level = 0.9)

committee.binary.table <- cbind(summary.committee_binary.95[[1]], summary.committee_binary.95[[2]], summary.committee_binary.9[[2]])
committee.binary <- xtable(committee.binary.table)
print(committee.binary, file = "butler_committee_binary_table4_topleft.txt")


# Weighted committee with copartisanship
summary.committee_weighted.95 <- BFP.results.summary(parameters, copartisan_committee_weighted_results, level = 0.95)
summary.committee_weighted.9 <- BFP.results.summary(parameters, copartisan_committee_weighted_results, level = 0.9)

committee.weighted.table <- cbind(summary.committee_weighted.95[[1]], summary.committee_weighted.95[[2]], summary.committee_weighted.9[[2]])
committee.weighted <- xtable(committee.weighted.table)
print(committee.weighted, file = "butler_committee_weighted_table4_topright.txt")



