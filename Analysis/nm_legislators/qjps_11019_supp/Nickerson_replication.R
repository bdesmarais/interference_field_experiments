rm(list=ls())

library(foreign) ## For importing STATA data tables
install.packages("scales")
library(scales) ## For rounding up numbers as required
install.packages("Zelig")
library(Zelig) ## For clustered standard errors
install.packages("xtable")
library(xtable) ## For creating tables that can be used in LaTex


data <- read.dta("nm.replication.dta", convert.underscore=TRUE)

sort(names(data))

data$dem.vs <- data$kerry / (data$bush + data$kerry)
data$v.richardson <- data$richard / (data$richard + data$dendahl)
data$spending.pref <- (1* data$q2full) + (0.3* data$q2reduced) + (0* data$q2none)
spmed <- as.numeric(quantile(data$spending.pref, 0.5))
data$prior.health <- data$healthsolutions == "Y"


########################################################
######## Creating table 1: randomization checks ########
########################################################

table1 <- matrix(NA, 8, 3)
table1 <- as.data.frame(table1)
colnames(table1) <- c("Treatment", "Control", "p-value")
rownames(table1) <- c("Republican",
               "Constituent support for spending",
               "Constituent support for health care",
               "Bush vote-share 04",
               "Member vote-share 06",
               "Running for re-election",
               "Running unopposed",
               "Supported prior health care bill")

data.treat <- data[data$treatment==1,]
data.control <- data[data$treatment==0,]
vars <- cbind(data$rep, data$q2full, data$q1favor, data$pres.2party, data$self.2party, data$runningforreelection, data$unopposed, data$prior.health)
vars.treat <- cbind(data.treat$rep, data.treat$q2full, data.treat$q1favor, data.treat$pres.2party, data.treat$self.2party, data.treat$runningforreelection, data.treat$unopposed, data.treat$prior.health)
vars.control <- cbind(data.control$rep, data.control$q2full, data.control$q1favor, data.control$pres.2party, data.control$self.2party, data.control$runningforreelection, data.control$unopposed, data.control$prior.health)


#### First 2 columns
table.1.means <- function(var){
  percent(mean(var, na.rm=TRUE))
}

table1[,1] <- apply(vars.treat, 2, FUN=table.1.means)
table1[,2] <- apply(vars.control, 2, FUN=table.1.means)


#### P-values in the table
pvals <- function(var){
  value <- summary(glm(data$treatment ~ var, family="binomial"))$coef[2,4]
  format(value, digits=2)
}

table1[,3] <- apply(vars, 2, FUN=pvals)
table1
sink("table1.for.LaTeX.txt", append=FALSE, split=TRUE)
toLatex.xtable(xtable(table1, caption="Table 1: Randomization checks"))
sink()
file.show("table1.for.LaTeX.txt")


#########################################################################
######## Creating table 2: Results with full spending preference ########
#########################################################################


q2med <- as.numeric(quantile(data$q2full, 0.5))
q2sd <- as.numeric(sd(data$q2full))
q225 <- as.numeric(quantile(data$q2full, 0.25))


data$full.co <- 0
for (i in 1:length(data$full.co)){
  if (data$q2full[i] <= q2med){
    data$full.co[i] <- 1
  }
}

data$flt <- data$treatment * data$full.co

## This is the model with treatment, low support for spending and their interaction
probit.no.control <- zelig(data$sb24 ~ data$treatment + data$full.co + data$flt,
                           data=data, model='probit', robust=T, cluster="data$match_category")
no.control.output <- as.table(summary(probit.no.control)$coef[,c(1:2)])
rownames(no.control.output) <- c("Constant", "Treatment", "Low support for spending", "Low support*Treatment")
colnames(no.control.output) <- c("Coefficients", "SEs")
no.control.output
format(no.control.output, digits=2)

sink("table2A.for.LaTeX.txt", append=FALSE, split=TRUE)
toLatex.xtable(xtable(no.control.output, caption="Table 2A: Regression results without controls"))
sink()
file.show("table2A.for.LaTeX.txt")

#write.table(no.control.output, file="T2_P1.csv", row.names=T, col.names=T, sep=",")


## This is the model with treatment, low support for spending and their interaction
## with controls for Partisanship and Vote Share
probit.with.control <- zelig(data$sb24 ~ data$treatment + data$full.co + data$flt + data$rep + data$dem.vs,
                           data=data, model='probit', robust=T, cluster="data$match_category")
summary(probit.with.control)
with.control.output <- as.table(summary(probit.with.control)$coef[,c(1:2)])
rownames(with.control.output) <- c("Constant", "Treatment",
                                 "Low support for spending", "Low support*Treatment",
                                 "Republican legislator", "2004 Democratic two-party presidential vote share")
colnames(with.control.output) <- c("Coefficients", "SEs")
with.control.output
format(with.control.output, digits=2)

sink("table2B.for.LaTeX.txt", append=FALSE, split=TRUE)
toLatex.xtable(xtable(with.control.output, caption="Regression results with controls"))
sink()
file.show("table2B.for.LaTeX.txt")


#write.table(with.control.output, file="T2_P2.csv", row.names=T, col.names=T, sep=",")


