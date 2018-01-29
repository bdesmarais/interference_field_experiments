# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 2 Analysis

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

# Bring in Data
load("Experiment_2_Permutation_Matrix.rdata")
load("Experiment_2_Past_Elections.rdata")
load("Experiment_2_Results.rdata")

# Clean and prepare data --------------------------------------------------

set.seed(12345)
chosen.rand <- sample(1:10000, 1)

albany <- within(albany,{
  restricted.prob.3 <- rowMeans(restricted.condition.block==3, na.rm=TRUE)
  restricted.prob.2 <- rowMeans(restricted.condition.block==2, na.rm=TRUE)
  restricted.prob.1 <- rowMeans(restricted.condition.block==1, na.rm=TRUE)
  prob.Z <- rowMeans(restricted.Z.block==1)
  Z.obs <- restricted.Z.block[,chosen.rand]
  condition.obs <- restricted.condition.block[,chosen.rand]
})

albany <- merge(albany, officialresults, by.x = "id", by.y = "ed")

albany <- within(albany,{
  kathyvoteshare <- kathy13/(kathy13+corey13+write13)
  kathyvotemargin <- kathy13 - (corey13+write13)
  turnout13 <- kathy13+corey13+write13
  
  share <- kathyvoteshare
  margin <- kathyvotemargin
  turnout <- turnout13
  
  jenningsvoteshare09 <- jennings09/(jennings09+ ellis09)
  jenningsvotemargin09 <- jennings09 - (ellis09)
  turnout09 <- jennings09+ ellis09
  
  share_cov <- jenningsvoteshare09
  margin_cov <- jenningsvotemargin09
  turnout_cov <- turnout09
  
  jenningsvoteshare05 <- jennings05/(jennings05+goodbee05)
  jenningsvotemargin05 <- jennings05 - (goodbee05)
  turnout05 <- jennings05+goodbee05
  
  
  weights <- rep(NA, length(id))
  weights[condition.obs==3] <- 1/restricted.prob.3[condition.obs==3]
  weights[condition.obs==2] <- 1/restricted.prob.2[condition.obs==2]
  weights[condition.obs==1] <- 1/restricted.prob.1[condition.obs==1]
  
  weights.Z <- rep(NA, length(id))
  weights.Z[Z.obs==1] <- 1/prob.Z[Z.obs==1]
  weights.Z[Z.obs==0] <- 1/(1-prob.Z[Z.obs==0])
  
  in_direct_exp <- as.numeric(condition.obs %in% c(1,3) & restricted.prob.3<1 & restricted.prob.3>0 & restricted.prob.1<1 & restricted.prob.1>0)
  in_indirect_exp <- as.numeric(condition.obs %in% c(1,2) & restricted.prob.2<1 & restricted.prob.2>0 & restricted.prob.1<1 & restricted.prob.1>0)
  condition_factor <- factor(condition.obs, levels=1:3, labels = c("Control", "Adjacent", "Treated"))
  condition_factor <- factor(condition_factor, levels=c("Control", "Treated", "Adjacent"))
  
  in_whole_exp <- as.numeric(restricted.prob.3<1 & restricted.prob.3>0 & 
                               restricted.prob.2<1 & restricted.prob.2>0 & 
                               restricted.prob.1<1 & restricted.prob.1>0)
  include <- in_whole_exp
  study <- "Experiment 2"
  
})

# Save the output to exp_2.rdata
# exp_2 <- albany
# cond_mat_2 <- restricted.condition.block
# save(exp_2,cond_mat_2, file = "exp_2.rdata")

## replication table two row two
table(albany$condition_factor, albany$include)

# Conduct Analysis --------------------------------------------------------

voteshare <- lm(share ~ condition_factor, weights=weights, data=subset(albany, in_whole_exp==1))
margin <- lm(margin ~ condition_factor, weights=weights, data=subset(albany, in_whole_exp==1))
turnout <- lm(turnout ~ condition_factor, weights=weights, data=subset(albany, in_whole_exp==1))

voteshare_cov <- lm(share ~ condition_factor + registereddems + jenningsvoteshare09 + jenningsvoteshare05, weights=weights, data=subset(albany, in_whole_exp==1))
margin_cov <- lm(margin ~ condition_factor + registereddems + jenningsvotemargin09 + jenningsvotemargin05, weights=weights, data=subset(albany, in_whole_exp==1))
turnout_cov <- lm(turnout ~ condition_factor + registereddems + turnout09 + turnout05, weights=weights, data=subset(albany, in_whole_exp==1))

# Export Tables -----------------------------------------------------------

## replication table four
# sink("experiment_2.tex")
stargazer(voteshare, voteshare_cov,
          se = makerobustseslist(list(voteshare, voteshare_cov)),
          title="Impact of Lawn Signs on Vote Share (Experiment 2)",
          label = "tab: exp_2",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2"),
          dep.var.labels= c("Vote Share"),
          covariate.labels=c("Assigned Lawn Signs (n=15)", "Adjacent to Lawn Signs (n=41)", NA),
          omit="(registereddems)|(05)|(09)", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{10cm}{Covariates: Registered Democrats and Mayoral Vote Margin '05 and '09.}"),
          notes.append=FALSE)
# sink()

## replication table c.2 appendix
# sink("experiment_2_appendix.tex")
stargazer(margin, margin_cov, turnout, turnout_cov,
          se = makerobustseslist(list(margin, margin_cov, turnout, turnout_cov)),
          title="Impact of Lawn Signs on Margin and Turnout (Experiment 2)",
          label = "tab: exp_2_app",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2", "Model 1", "Model 2"),
          dep.var.labels= c("Vote Margin", "Turnout"),
          covariate.labels=c("Assigned Lawn Signs (n=15)", "Adjacent to Lawn Signs (n=41)", NA),
          omit="(registereddems)|(05)|(09)", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{13cm}{Column 2 covariates: Registered Democrats and Mayoral Vote Share '05 and '09. Column 4 covariates: Registered Democrats and Mayoral Turnout '05 and '09.}"),
          notes.append=FALSE)
# sink()

# Randomization Inference -------------------------------------------------

f_margin <-
  get_f(formula_F = "margin ~ condition_factor + registereddems + jenningsvotemargin09 + jenningsvotemargin05",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = subset(albany, in_whole_exp==1))

f_share <-
  get_f(formula_F = "share ~ condition_factor + registereddems + jenningsvoteshare09 + jenningsvoteshare05",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = subset(albany, in_whole_exp==1))

f_turnout <-
  get_f(formula_F = "turnout ~ condition_factor + registereddems + turnout09 + turnout05",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = subset(albany, in_whole_exp==1))

sims<- 10000
f_margin_sims <- rep(NA, sims)
f_share_sims <- rep(NA, sims)
f_turnout_sims <- rep(NA, sims)

set.seed(12345)
for(i in 1:sims){
  condition.sim <- restricted.condition.block[,i]
  weights.sim <- rep(NA, length(condition.sim))
  weights.sim[condition.sim==3] <- 1/(albany$restricted.prob.3[condition.sim==3])
  weights.sim[condition.sim==2] <- 1/(albany$restricted.prob.2[condition.sim==2])
  weights.sim[condition.sim==1] <- 1/(albany$restricted.prob.1[condition.sim==1])
  
  condition_sim_factor <- factor(condition.sim, levels = 1:3, labels = c("Control", "Adjacent", "Treated"))
  condition_sim_factor <- factor(condition_sim_factor, levels=c("Control", "Treated", "Adjacent"))
  
  frame.sim <- data.frame(albany, condition_sim_factor, weights.sim)
  
  f_margin_sims[i] <-
    get_f(formula_F = "margin ~ condition_sim_factor + registereddems + jenningsvotemargin09 + jenningsvotemargin05",
          treat_var = "condition_sim_factor", 
          weights = "weights.sim", 
          data = subset(frame.sim, in_whole_exp==1))
  
  f_share_sims[i] <-
    get_f(formula_F = "share ~ condition_sim_factor + registereddems + jenningsvoteshare09 + jenningsvoteshare05",
          treat_var = "condition_sim_factor", 
          weights = "weights.sim", 
          data = subset(frame.sim, in_whole_exp==1))
  
  f_turnout_sims[i] <-
    get_f(formula_F = "turnout ~ condition_sim_factor + registereddems + turnout09 + turnout05",
          treat_var = "condition_sim_factor", 
          weights = "weights.sim", 
          data = subset(frame.sim, in_whole_exp==1))
  
  print(i)
}


pval_share<-mean(f_share_sims >= f_share)
pval_margin <-mean(f_margin_sims >= f_margin)
pval_turnout <-mean(f_turnout_sims >= f_turnout)

ps <- t(matrix(c(pval_margin, pval_share, pval_turnout)))
colnames(ps) <- c("Margin", "Share",  "Turnout")
rownames(ps) <- "p(f_sims >= f_obs)"

# sink(file = "experiment_2_ps.tex")
## replication randomization pvalue referred to in text, section 5.2
xtable(ps, caption = "Experiment 2 Randomization Inference p-values")
# sink()
