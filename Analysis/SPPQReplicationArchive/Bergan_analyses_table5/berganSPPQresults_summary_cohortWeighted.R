# read in results

all_parameters <- NULL
copartisan_cohort_weighted_results <- NULL


load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity1.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity2.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity3.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity4.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity5.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity6.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity7.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity8.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity9.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity10.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity11.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity12.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity13.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity14.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity15.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_copartisan_cohort_chamber_similarity16.RData")
copartisan_cohort_weighted_results <- c(copartisan_cohort_weighted_results,do.call('rbind',BFP.results)[,1])
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

summary.copartisan_cohort_weighted <- BFP.results.summary(parameters,copartisan_cohort_weighted_results)
summary.copartisan_cohort_weighted.9 <- BFP.results.summary(parameters,copartisan_cohort_weighted_results,level=0.9)


library(xtable)

cohort.table <- cbind(summary.copartisan_cohort_weighted[[1]],summary.copartisan_cohort_weighted[[2]],summary.copartisan_cohort_weighted.9[[2]])
cohort.weighted <- xtable(cohort.table)
print(cohort.weighted, file = "table5_bottomright_bergan_cohort_weighted.txt")
