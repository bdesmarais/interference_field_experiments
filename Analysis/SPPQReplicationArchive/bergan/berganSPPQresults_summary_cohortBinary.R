# read in results

all_parameters <- NULL
copartisan_cohort_binary_results <- NULL


load("BerganSPPQRRresults_copartisan_cohort_chamber_binary1.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary2.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary3.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary4.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary5.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary6.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary7.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary8.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary9.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary10.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary11.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary12.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary13.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary14.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary15.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_binary16.RData")
copartisan_cohort_binary_results <- c(copartisan_cohort_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

parameters <- all_parameters


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

summary.copartisan_cohort_binary <- BFP.results.summary(parameters,copartisan_cohort_binary_results)
summary.copartisan_cohort_binary.9 <- BFP.results.summary(parameters,copartisan_cohort_binary_results,level=0.9)


library(xtable)

cohort.table <- cbind(summary.copartisan_cohort_binary[[1]],summary.copartisan_cohort_binary[[2]],summary.copartisan_cohort_binary.9[[2]])
xtable(cohort.table)
