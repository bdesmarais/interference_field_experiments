# read in results

all_parameters <- NULL
cospon_binary_results <- NULL


load("BerganSPPQRRresults_cospon_binary1.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary2.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary3.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary4.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary5.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary6.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary7.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary8.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary9.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary10.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary11.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary12.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary13.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary14.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary15.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
all_parameters <- rbind(all_parameters,parameters)

load("BerganSPPQRRresults_cospon_binary16.RData")
cospon_binary_results <- c(cospon_binary_results,do.call('rbind',BFP.results)[,1])
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

summary.cospon_binary <- BFP.results.summary(parameters,cospon_binary_results)
summary.cospon_binary.9 <- BFP.results.summary(parameters,cospon_binary_results,level=0.9)


library(xtable)

cospon.table <- cbind(summary.cospon_binary[[1]],summary.cospon_binary[[2]],summary.cospon_binary.9[[2]])
xtable(cospon.table)
