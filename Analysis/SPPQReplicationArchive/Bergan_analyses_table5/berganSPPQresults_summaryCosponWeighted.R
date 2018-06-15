# read in results
load("BerganSPPQRRresults_cospon_weighted_newmodel.RData")
copartisan_cosponsorship_weighted_results <- do.call('rbind',BFP.results)[,1]


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

summary.copartisan_cospon_weighted <- BFP.results.summary(parameters,copartisan_cosponsorship_weighted_results)

summary.copartisan_cospon_weighted.9 <- BFP.results.summary(parameters,copartisan_cosponsorship_weighted_results,level=0.9)

library(xtable)

cospon.table <- cbind(summary.copartisan_cospon_weighted[[1]],summary.copartisan_cospon_weighted[[2]],summary.copartisan_cospon_weighted.9[[2]])
cospon.weighted <- xtable(cospon.table)
print(cospon.weighted, file = "table5_topright_bergan_cospon_weighted.txt")

