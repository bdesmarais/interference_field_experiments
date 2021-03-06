# read in results
load("CoppockSPPQRRresults_copartisan_committee_binary.RData")
copartisan_committee_binary_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_copartisan_committee_weighted.RData")
copartisan_committee_weighted_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_threshold60_unscaled.RData")
threshold60_unscaled_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_threshold60_scaled.RData")
threshold60_scaled_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_threshold80_scaled.RData")
threshold80_scaled_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_threshold80_unscaled.RData")
threshold80_unscaled_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_copartisan_cohort_binary.RData")
copartisan_cohort_binary_results <- do.call('rbind',BFP.results)[,1]

load("CoppockSPPQRRresults_copartisan_cohort_weighted.RData")
copartisan_cohort_weighted_results <- do.call('rbind',BFP.results)[,1]



BFP.results.summary <- function(parameters,p.values,level=0.95){
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

summary.copartisan_committee_binary <- BFP.results.summary(parameters,copartisan_committee_binary_results)
summary.copartisan_committee_weighted <- BFP.results.summary(parameters,copartisan_committee_weighted_results)
summary.copartisan_cohort_binary <- BFP.results.summary(parameters, copartisan_cohort_binary_results)
summary.copartisan_cohort_weighted <- BFP.results.summary(parameters, copartisan_cohort_weighted_results)
summary.copartisan_committee_binary.9 <- BFP.results.summary(parameters,copartisan_committee_binary_results,level=0.9)
summary.copartisan_committee_weighted.9 <- BFP.results.summary(parameters,copartisan_committee_weighted_results,level=0.9)
summary.copartisan_cohort_binary.9 <- BFP.results.summary(parameters, copartisan_cohort_binary_results,level=0.9)
summary.copartisan_cohort_weighted.9 <- BFP.results.summary(parameters, copartisan_cohort_weighted_results,level=0.9)

library(xtable)

committee.table <- cbind(summary.copartisan_committee_binary[[1]],summary.copartisan_committee_binary[[2]],summary.copartisan_committee_binary.9[[2]],summary.copartisan_committee_weighted[[1]],summary.copartisan_committee_weighted[[2]],summary.copartisan_committee_weighted.9[[2]])
cohort.table <- cbind(summary.copartisan_cohort_binary[[1]],summary.copartisan_cohort_binary[[2]],summary.copartisan_cohort_binary.9[[2]],summary.copartisan_cohort_weighted[[1]],summary.copartisan_cohort_weighted[[2]],summary.copartisan_cohort_weighted.9[[2]])
xtable(committee.table)
