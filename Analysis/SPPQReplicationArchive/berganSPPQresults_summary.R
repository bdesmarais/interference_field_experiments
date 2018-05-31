####
#### Creating tables using results
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


## Read in results of all four Bergan analyses
# Binary cohort with copartisanship
load("BerganSPPQRRresults_copartisan_cohort_binary.RData")
copartisan_cohort_binary_results <- do.call('rbind',BFP.results)[,1]

# Weighted cohort with copartisanship
load("BerganSPPQRRresults_copartisan_cohort_similarity.RData")
copartisan_cohort_weighted_results <- do.call('rbind',BFP.results)[,1]

# Binary cosponsorship with copartisanship
load("BerganSPPQRRresults_cospon_binary.RData")
copartisan_cosponsorship_binary_results <- do.call('rbind',BFP.results)[,1]

# Weighted cosponsorship with copartisanship
load("BerganSPPQRRresults_cospon_weighted.RData")
copartisan_cosponsorship_weighted_results <- do.call('rbind',BFP.results)[,1]


## Generate the point estimates and 90 & 95% CIs for each analysis
# Binary cohort with copartisanship
summary.cohort_binary.95 <- BFP.results.summary(parameters, copartisan_cohort_binary_results, level = 0.95)
summary.cohort_binary.9 <- BFP.results.summary(parameters, copartisan_cohort_binary_results,level = 0.9)

cohort.binary.table <- cbind(summary.cohort_binary.95[[1]], summary.cohort_binary.95[[2]], summary.cohort_binary.9[[2]])
xtable(cohort.binary.table)


# Weighted cohort with copartisanship
summary.cohort_weighted.95 <- BFP.results.summary(parameters, copartisan_cohort_weighted_results, level = 0.95)
summary.cohort_weighted.9 <- BFP.results.summary(parameters, copartisan_cohort_weighted_results,level = 0.9)

cohort.weighted.table <- cbind(summary.cohort_weighted.95[[1]], summary.cohort_weighted.95[[2]], summary.cohort_weighted.9[[2]])
xtable(cohort.weighted.table)


# Binary cosponsorship with copartisanship
summary.cospon_binary.95 <- BFP.results.summary(parameters, copartisan_cosponsorship_binary_results, level = 0.95)
summary.cospon_binary.9 <- BFP.results.summary(parameters, copartisan_cosponsorship_binary_results,level = 0.9)

cospon.binary.table <- cbind(summary.cospon_binary.95[[1]], summary.cospon_binary.95[[2]], summary.cospon_binary.9[[2]])
xtable(cospon.binary.table)


# Weighted cosponsorship with copartisanship
summary.cospon_weighted.95 <- BFP.results.summary(parameters, copartisan_cosponsorship_weighted_results, level = 0.95)
summary.cospon_weighted.9 <- BFP.results.summary(parameters, copartisan_cosponsorship_weighted_results, level = 0.9)

cospon.weighted.table <- cbind(summary.cospon_weighted.95[[1]], summary.cospon_weighted.95[[2]], summary.cospon_weighted.9[[2]])
xtable(cospon.weighted.table)



