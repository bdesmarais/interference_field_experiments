# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 1 Analysis

# Uncomment to set your working directory
# setwd("")

# Uncomment to install required packages
# install.packages(c("stargazer", "xtable", "sandwich", "lmtest"))

rm(list=ls())

# Load Packages, Data, and Functions --------------------------------------

library(stargazer)
library(xtable)

# Helper functions
source("Lawn_Signs_Source.R")

# Read in Adjacency Data
load("Experiment_1_Adjacencies.rdata")

# Read in probabilites of assignment and 10K potential randomizations
load("Experiment_1_Permutation_Matrix.rdata")

# Read in voting results
load("Experiment_1_Results.rdata")

# Clean and prepare data --------------------------------------------------

exp_1 <- within(exp_1, {
                      dvote12 <- c12demvotes/(c12demvotes+ c12repvotes+c12indvotes)
                      dvote10 <- c10demvotes/(c10demvotes+ c10repvotes+c10indvotes)
                      dvote08 <- c08demvotes/(c08demvotes+ c08repvotes+c08indvotes+c08convotes)
                      dvote06 <- c06demvotes/(c06demvotes+ c06repvotes+c06indvotes+c06convotes)
                      pvote08 <- (p08demvotes+ p08wfpvotes)/(p08demvotes+ p08repvotes+p08indvotes+p08convotes+ p08wfpvotes)
                      pvote12 <- (p12demvotes+ p12wfpvotes)/ (p12demvotes +p12repvotes +p12convotes +p12wfpvotes +p12greenvotes +p12socvotes +  p12lbnvotes + p12constvotes + p12othervotes)
                      
                      #vote margin
                      nvote12 <- c12demvotes - (c12repvotes+c12indvotes)
                      npvote12 <- (p12demvotes+ p12wfpvotes) - (p12repvotes +p12convotes) ## mitt wasn't a member of the independence party
                      nvote10 <- c10demvotes - (c10repvotes+c10indvotes)
                      nvote08 <- c08demvotes - (c08repvotes+c08indvotes+c08convotes)
                      nvote06 <- c06demvotes - (c06repvotes+c06indvotes+c06convotes)
                      npvote08 <- p08demvotes+ p08wfpvotes - (p08repvotes+p08indvotes+p08convotes)
                     
                     #total votes cast
                      totalcvotes12 <- c12demvotes + c12repvotes+ c12indvotes+ c12othrvotes+ c12voidvotes+ c12blankvotes 
                      totalcvotes10 <- c10demvotes + c10repvotes+ c10indvotes+ c10othrvotes+ c10voidvotes+ c10blankvotes 
                      totalcvotes08 <- c08demvotes + c08repvotes+ c08indvotes+ c08othrvotes+ c08voidvotes+ c08blankvotes 
                      totalcvotes06 <- c06demvotes + c06repvotes+ c06indvotes+ c06othrvotes+ c06voidvotes+ c06blankvotes
                      totalpvotes08 <- p08demvotes + p08wfpvotes+ p08repvotes+ p08indvotes+ p08convotes
                      
                      share <- dvote12
                      margin <- nvote12
                      turnout <- totalcvotes12
                      
                      share_cov <- dvote10
                      margin_cov <- nvote10
                      turnout_cov <- totalcvotes10
})

# this step aggregates the data into the merged electoral districts

lawnresults <- aggregate(exp_1[,c(3, 143:165)], by=list(exp_1$fid), sum)  
names(lawnresults)[1] <- "fid"
lawnresults <- lawnresults[order(lawnresults$fid),]

untreatables <- c("11", "79", "80", "82", "85")  ### These are all the untreated except 70, which we don't think is parallel.

lawnresults <- within(lawnresults,{
                      Z.obs <- selected
                      Z.obs[fid %in% untreatables] <- 0
                      treatable <- rep(1, 93)
                      treatable[fid %in% untreatables] <- 0
                    
                     ## Determining observed treatment condition (where "Z.obs" is the same as treated)
                      numberofspillovers <- A.93 %*% Z.obs
                      condition <- rep(NA, 93)
                      condition[Z.obs==1 & numberofspillovers>0] <- 1
                      condition[Z.obs==1 & numberofspillovers==0] <- 2
                      condition[Z.obs==0 & numberofspillovers>0] <- 3
                      condition[Z.obs==0 & numberofspillovers==0] <- 4
                      
                    ### Giving observations the weights for the conditions they are actually in.
                      weights <- rep(NA, 93)
                             for (i in 1:length(weights)){
                             if(condition[i]==1){weights[i] <- 1/probs.1.1[i]}
                             if(condition[i]==2){weights[i] <- 1/probs.1.0[i]}
                             if(condition[i]==3){weights[i] <- 1/probs.0.1[i]}
                             if(condition[i]==4){weights[i] <- 1/probs.0.0[i]}
                             }
                      condition_factor <- factor(condition, levels=2:4, labels = c("Treated", "Adjacent", "Control"))
                      condition_factor <- relevel(condition_factor, ref="Control")
                      probs.0.0 <- probs.0.0
                      probs.0.1 <- probs.0.1
                      probs.1.0 <- probs.1.0
                      probs.1.1 <- probs.1.1
                      include <- treatable
                      study <- "Experiment 1"
})



# Save the output to exp_1.rdata
# exp_1 <- lawnresults
# cond_mat_1 <- condition.block
# save(exp_1,cond_mat_1, file = "exp_1.rdata")

## replication table two
table(lawnresults$condition_factor, lawnresults$treatable)

# Conduct Analysis --------------------------------------------------------

# Congressional Vote Margin
nvote12 <- lm(nvote12 ~ condition_factor, weights=weights, data=subset(lawnresults, treatable==1))
cov.nvote12 <- lm(nvote12 ~ condition_factor + nvote10 + nvote08 + nvote06 + npvote08, weights=weights, 
                  data=subset(lawnresults, treatable==1))

# Congressional Vote Share
dvote12 <- lm(dvote12 ~ condition_factor, weights=weights, data=subset(lawnresults, treatable==1))
cov.dvote12 <- lm(dvote12 ~ condition_factor + dvote10 + dvote08 + dvote06 + pvote08, weights=weights, 
                  data=subset(lawnresults, treatable==1))

# Congressional Votes Cast
total12 <- lm(totalcvotes12 ~ condition_factor, weights=weights, data=subset(lawnresults, treatable==1))
cov.total12 <- lm(totalcvotes12 ~ condition_factor + totalcvotes10 + totalcvotes08 + totalcvotes06 + totalpvotes08,
                  weights=weights, data=subset(lawnresults, treatable==1))

# Export Tables -----------------------------------------------------------

library(stargazer)
## replication table three
# sink("experiment_1.tex")
stargazer(dvote12, cov.dvote12,
          se = makerobustseslist(list(dvote12, cov.dvote12)),
          title="Impact of Lawn Signs on Vote Share (Experiment 1)",
          label = "tab: exp_1",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2"),
          dep.var.labels= c("Vote Share"),
          covariate.labels=c("Assigned Lawn Signs (n=23)", "Adjacent to Lawn Signs (n=49)", NA),
          omit="vote", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{10cm}{Covariates: Congressional Vote Margin '06, '08, '10 and Presidential Vote Margin '08.}"), 
          notes.append=FALSE)
# sink()

## replication table C.1 appendix
# sink("experiment_1_appendix.tex")
stargazer(nvote12, cov.nvote12, total12, cov.total12,
          se = makerobustseslist(list(nvote12, cov.nvote12, total12, cov.total12)),
          title="Impact of Lawn Signs on Margin and Turnout (Experiment 1)",
          label = "tab: exp_1_app",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2", "Model 1", "Model 2"),
          dep.var.labels= c("Vote Margin", "Turnout"),
          covariate.labels=c("Assigned Lawn Signs (n=23)", "Adjacent to Lawn Signs (n=49)", NA),
          omit="vote", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{12.5cm}{Column 2 covariates: Congressional Vote Share '06, '08, '10 and Presidential Vote Share '08. Column 4 covariates: Congressional Turnout '06, '08, '10 and Presidential Turnout '08}"),
          notes.append=FALSE)
# sink()

# Randomization Inference -------------------------------------------------

f_margin <-
  get_f(formula_F = "nvote12 ~ condition_factor + nvote10 + nvote08 + nvote06 + npvote08",
      treat_var = "condition_factor", 
      weights = "weights", 
      data = subset(lawnresults, treatable==1))

f_share <-
  get_f(formula_F = "dvote12 ~ condition_factor + dvote10 + dvote08 + dvote06 + pvote08",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = subset(lawnresults, treatable==1))

f_turnout <-
  get_f(formula_F = "totalcvotes12 ~ condition_factor + totalcvotes10 + totalcvotes08 + totalcvotes06 + totalpvotes08",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = subset(lawnresults, treatable==1))

sims <- 10000
f_margin_sims <- rep(NA, sims)
f_share_sims <- rep(NA, sims)
f_turnout_sims <- rep(NA, sims)

for (i in 1:sims){
  #print(i) 
  condition.sims <- condition.block[,i]
  weights.sims <- rep(NA, 93)
  
  for (j in 1:93){
    if(condition.sims[j]==1){weights.sims[j] <- 1/probs.1.1[j]}
    if(condition.sims[j]==2){weights.sims[j] <- 1/probs.1.0[j]}
    if(condition.sims[j]==3){weights.sims[j] <- 1/probs.0.1[j]}
    if(condition.sims[j]==4){weights.sims[j] <- 1/probs.0.0[j]}
  }
  
  condition_factor_sims <- factor(condition.sims, levels=2:4, labels = c("Treated", "Adjacent", "Control"))
  condition_factor_sims <- relevel(condition_factor_sims, ref="Control")
  lawnresults.sims <- data.frame(lawnresults, condition_factor_sims, weights.sims)
  f_margin_sims[i] <-
    get_f(formula_F = "nvote12 ~ condition_factor_sims + nvote10 + nvote08 + nvote06 + npvote08",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = subset(lawnresults.sims, treatable==1))
  
  f_share_sims[i] <-
    get_f(formula_F = "dvote12 ~ condition_factor_sims + dvote10 + dvote08 + dvote06 + pvote08",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = subset(lawnresults.sims, treatable==1))
  
  f_turnout_sims[i] <-
    get_f(formula_F = "totalcvotes12 ~ condition_factor_sims + totalcvotes10 + totalcvotes08 + totalcvotes06 + totalpvotes08",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = subset(lawnresults.sims, treatable==1))
}

pval_share<-mean(f_share_sims >= f_share)
pval_margin <-mean(f_margin_sims >= f_margin)
pval_turnout <-mean(f_turnout_sims >= f_turnout)

ps <- t(matrix(c(pval_margin, pval_share, pval_turnout)))
colnames(ps) <- c("Margin", "Share",  "Turnout")
rownames(ps) <- "p(f_sims >= f_obs)"

# sink(file ="experiment_1_ps.tex")
## replication randomization pvalue referred to in text, section 5.1
xtable(ps, caption = "Experiment 1 Randomization Inference p-values")
# sink()
