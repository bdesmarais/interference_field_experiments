# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 3 Analysis

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

# Bring in data
load("Experiment_3_Data_with_Random_Assignment.rdata")
load("Experiment_3_Results.rdata")

# Clean and prepare data --------------------------------------------------

results <- merge(x=results, y=results2013, by.x="precinct_id_2012", by.y = "precinct_id")

results <- within(results,{
  writein[is.na(writein)] <- 0
  
  cuccinelli_margin_2013 <- Cuccinelli - (McAuliffe + Sarvis + writein)
  cuccinelli_share_2013 <- Cuccinelli/ (Cuccinelli + McAuliffe + Sarvis + writein)
  cuccinelli_2margin_2013 <- Cuccinelli - (McAuliffe)
  cuccinelli_2share_2013 <- Cuccinelli/ (Cuccinelli + McAuliffe)
  
  McAuliffe_margin_2013 <- McAuliffe - (Cuccinelli + Sarvis + writein)
  McAuliffe_share_2013 <- McAuliffe/ (Cuccinelli + McAuliffe + Sarvis + writein)
  McAuliffe_2margin_2013 <- McAuliffe - Cuccinelli
  McAuliffe_2share_2013 <- McAuliffe/ (Cuccinelli + McAuliffe)
  
  turnout_2013 <-  Cuccinelli + McAuliffe + Sarvis + writein
    
  obama_margin_2012 <- Barack.Obama..D..2012.Vote.Totals - (Mitt.Romney..R..2012.Vote.Totals + OTHER.2012.Vote.Totals)
  obama_share_2012 <- Barack.Obama..D..2012.Vote.Totals/ (Barack.Obama..D..2012.Vote.Totals + OTHER.2012.Vote.Totals + Mitt.Romney..R..2012.Vote.Totals)
  obama_2margin_2012 <- Barack.Obama..D..2012.Vote.Totals - (Mitt.Romney..R..2012.Vote.Totals)
  obama_2share_2012 <- Barack.Obama..D..2012.Vote.Totals/ (Barack.Obama..D..2012.Vote.Totals + Mitt.Romney..R..2012.Vote.Totals)
  
  turnout_2012 <- Barack.Obama..D..2012.Vote.Totals + OTHER.2012.Vote.Totals + Mitt.Romney..R..2012.Vote.Totals
  
  shannon_margin_2009 <- Stephen.C..Shannon..D..Vote.Totals - (Ken.T..Cuccinelli.II..R..Vote.Totals   + Other.2009.Vote.Totals)
  shannon_share_2009 <- Stephen.C..Shannon..D..Vote.Totals / (Stephen.C..Shannon..D..Vote.Totals   + Other.2009.Vote.Totals + Ken.T..Cuccinelli.II..R..Vote.Totals)

  shannon_2margin_2009 <- Stephen.C..Shannon..D..Vote.Totals - (Ken.T..Cuccinelli.II..R..Vote.Totals)
  shannon_2share_2009 <- Stephen.C..Shannon..D..Vote.Totals / (Stephen.C..Shannon..D..Vote.Totals + Ken.T..Cuccinelli.II..R..Vote.Totals)
  
  turnout_2009 <- Stephen.C..Shannon..D..Vote.Totals   + Other.2009.Vote.Totals + Ken.T..Cuccinelli.II..R..Vote.Totals
  
  margin <- cuccinelli_margin_2013
  share <- cuccinelli_share_2013
  turnout <- turnout_2013
  
  margin_cov <- obama_margin_2012
  share_cov <- obama_share_2012
  turnout_cov <- turnout_2012
  
  
  probs.Z <- rowMeans(restricted.Z.block)
  weights.Z <- rep(NA, 131)
  weights.Z[Z.obs==1] <- 1/probs.Z[Z.obs==1]
  weights.Z[Z.obs==0] <- 1/(1- probs.Z[Z.obs==0])
  
  condition_factor <- factor(condition.obs, levels = 1:3, labels = c("Control", "Adjacent", "Treated"))
  condition_factor <- factor(condition_factor, levels=c("Control", "Treated", "Adjacent"))
  
  include <- 1
  study <- "Experiment 3"
})

# Save the output to exp_3.rdata
# exp_3 <- results
# cond_mat_3 <- restricted.condition.block
# save(exp_3, cond_mat_3, file = "exp_3.rdata")
## replication table two row three
table(results$condition_factor, results$include)

# Conduct Analysis --------------------------------------------------------

margin.fit <- lm(margin ~  condition_factor, weights=weights, data=results)
margin.fit.cov <- lm(margin ~  condition_factor + obama_margin_2012 + shannon_margin_2009 , weights=weights, data=results)

share.fit <- lm(share ~  condition_factor, weights=weights, data=results)
share.fit.cov <- lm(share ~  condition_factor + obama_share_2012 + shannon_share_2009 , weights=weights, data=results)

turnout.2013.fit <- lm(turnout ~  condition_factor, weights=weights, data=results)
turnout.2013.fit.cov <- lm(turnout ~  condition_factor + turnout_2012 + turnout_2009, weights=weights, data=results)

# Export Tables -----------------------------------------------------------

## replication table five
# sink("experiment_3.tex")
stargazer(share.fit, share.fit.cov,
          se = makerobustseslist(list(share.fit, share.fit.cov)),
          title="Impact of Lawn Signs on Vote Share (Experiment 3)",
          label = "tab: exp_3",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2"),
          dep.var.labels= c("Vote Share"),
          covariate.labels=c("Assigned Lawn Signs (n=30)", "Adjacent to Lawn Signs (n=76)", NA),
          omit=c("(2012)|(2009)"), 
          omit.labels = c("Covariate Adjustment"),
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{10cm}{Covariates: Gubernatorial Vote Margin '09 and Presidential Vote Margin '12.}"),
          notes.append=FALSE)
# sink()

## replication appendix table C.3
# sink("experiment_3_appendix.tex")
stargazer(margin.fit, margin.fit.cov, turnout.2013.fit, turnout.2013.fit.cov,
          se = makerobustseslist(list(margin.fit, margin.fit.cov, turnout.2013.fit, turnout.2013.fit.cov)),
          title="Impact of Lawn Signs on Margin and Turnout (Experiment 3)",
          label = "tab: exp_3_app",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2", "Model 1", "Model 2"),
          dep.var.labels= c("Vote Margin", "Turnout"),
          covariate.labels=c("Assigned Lawn Signs (n=30)", "Adjacent to Lawn Signs (n=76)", NA),
          omit=c("(2012)|(2009)"), 
          omit.labels = c("Covariate Adjustment"),
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{13cm}{Column 2 covariates: Gubernatorial Vote Share '09 and Presidential Vote Share '12. Column 4 covariates: Gubernatorial Vote Turnout '09 and Presidential Vote Turnout '12.}"),
          notes.append=FALSE)
# sink()

# Randomization Inference -------------------------------------------------

f_margin <-
  get_f(formula_F = "margin ~  condition_factor + obama_margin_2012 + shannon_margin_2009",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = results)

f_share <-
  get_f(formula_F = "share ~  condition_factor + obama_share_2012 + shannon_share_2009 ",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = results)

f_turnout <-
  get_f(formula_F = "turnout ~  condition_factor + turnout_2012 + turnout_2009",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = results)



sims<- 10000
f_margin_sims <- rep(NA, sims)
f_share_sims <- rep(NA, sims)
f_turnout_sims <- rep(NA, sims)

set.seed(12345)
for(i in 1:sims){
condition.sim <- restricted.condition.block[,i]
weights.sim <- rep(NA, length(condition.sim))
weights.sim[condition.sim==3] <- 1/(results$restricted.prob.3[condition.sim==3])
weights.sim[condition.sim==2] <- 1/(results$restricted.prob.2[condition.sim==2])
weights.sim[condition.sim==1] <- 1/(results$restricted.prob.1[condition.sim==1])

condition_sim_factor <- factor(condition.sim, levels = 1:3, labels = c("Control", "Adjacent", "Treated"))
condition_sim_factor <- factor(condition_sim_factor, levels=c("Control", "Treated", "Adjacent"))
frame.sim <- data.frame(results, condition_sim_factor, weights.sim)

f_margin_sims[i] <-
  get_f(formula_F = "margin ~  condition_sim_factor + obama_margin_2012 + shannon_margin_2009",
        treat_var = "condition_sim_factor", 
        weights = "weights.sim", 
        data = frame.sim)

f_share_sims[i] <-
  get_f(formula_F = "share ~  condition_sim_factor + obama_share_2012 + shannon_share_2009 ",
        treat_var = "condition_sim_factor", 
        weights = "weights.sim", 
        data = frame.sim)

f_turnout_sims[i] <-
  get_f(formula_F = "turnout ~  condition_sim_factor + turnout_2012 + turnout_2009",
        treat_var = "condition_sim_factor", 
        weights = "weights.sim", 
        data = frame.sim)

print(i)
}


pval_share<-mean(f_share_sims >= f_share)
pval_margin <-mean(f_margin_sims >= f_margin)
pval_turnout <-mean(f_turnout_sims >= f_turnout)

ps <- t(matrix(c(pval_margin, pval_share, pval_turnout)))
colnames(ps) <- c("Margin", "Share",  "Turnout")
rownames(ps) <- "p(f_sims >= f_obs)"

#sink(file ="experiment_3_ps.tex")
## replication randomization pvalue referred to in text, section 5.3
xtable(ps, caption = "Experiment 3 Randomization Inference p-values")
#sink()
