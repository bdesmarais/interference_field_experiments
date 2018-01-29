# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 4 Analysis

# Uncomment to set your working directory
# setwd("")

# Uncomment to install required packages
# install.packages(c("stargazer", "xtable", "sandwich", "lmtest", "dplyr"))

rm(list=ls())

# Load Packages, Data, and Functions --------------------------------------

library(stargazer)
library(xtable)
library(dplyr)

# Helper functions
source("Lawn_Signs_Source.R")

# Bring in data
load("Experiment_4_Data_with_Random_Assignment.rdata")
load("Experiment_4_Permutation_Matrix.rdata")
load("Experiment_4_Results.rdata")

# Clean and prepare data --------------------------------------------------

PA_41_df <- 
  left_join(PA_41_df,results ) %>%
  mutate(
        turnout = Cross + Difilippo + Eichelberger + Schin + Other,
        share = (Eichelberger + Schin)/turnout,
        in_study = res_prob_treated<1 & res_prob_treated>0 & 
          res_prob_adj<1 & res_prob_adj>0 & 
          res_prob_control<1 & res_prob_control>0
  ) %>%
  filter(in_study) 

PA_41_df <- 
  mutate(PA_41_df,  
    weights = 1/get_prob_obs(condition_obs),
    condition_obs = factor(condition_obs),
    condition_obs = relevel(condition_obs, ref = "control"),
    margin = Eichelberger - (Cross + Difilippo + Schin + Other),
    GOVT2010 = GOVRV2010 + GOVDV2010,
    USPT2008 = USPRV2008 + USPDV2008,
    GOVT2006 = GOVRV2006 + GOVDV2006,
    USPT2004 = USPRV2004 + USPDV2004,
    GOVT2002 = GOVRV2002 + GOVDV2002,
    USPT2000 = USPRV2000 + USPDV2000,
    GOVM2010 = GOVRV2010 - GOVDV2010,
    USPM2008 = USPRV2008 - USPDV2008,
    GOVM2006 = GOVRV2006 - GOVDV2006,
    USPM2004 = USPRV2004 - USPDV2004,
    GOVM2002 = GOVRV2002 - GOVDV2002,
    USPM2000 = USPRV2000 - USPDV2000,
    margin_cov = GOVM2010,
    share_cov = GOVP2010,
    turnout_cov = GOVT2010)

PA_41_df <- within(PA_41_df, {
  condition_factor <- factor(condition_obs, levels=c("control", "adjacent", "treated"), labels = c("Control", "Adjacent", "Treated"))
  condition_factor <- factor(condition_factor, levels=c("Control", "Treated", "Adjacent"))
  include <- 1
  study <- "Experiment 4"
})

# Save the output to exp_4.rdata
# exp_4 <- PA_41_df
# cond_mat_4 <- restricted_condition_mat
# save(exp_4,cond_mat_4, file = "exp_4.rdata")

## replication table two row four
table(PA_41_df$condition_factor, PA_41_df$include)

# Conduct Analysis --------------------------------------------------------

fit_margin <- lm(margin ~ condition_factor, weights = weights, data=PA_41_df)
fit_cov_margin <- lm(margin ~ condition_factor +
                       GOVM2010 + USPM2008 + GOVM2006 + USPM2004 + GOVM2002 + USPM2000,
                     weights = weights, data=PA_41_df)

fit_share <- lm(share ~ condition_factor, weights = weights, data=PA_41_df)
fit_cov_share <- lm(share ~ condition_factor +
                      GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000,
                    weights = weights, data=PA_41_df)
fit_turnout <- lm(turnout ~ condition_factor, weights = weights, data=PA_41_df)
fit_cov_turnout <- lm(turnout ~ condition_factor +
                        GOVT2010 + USPT2008 + GOVT2006 + USPT2004 + GOVT2002 + USPT2000,
                      weights = weights, data=PA_41_df)

# Export Tables -----------------------------------------------------------

# sink("experiment_4.tex")
## replication table six
stargazer(fit_share, fit_cov_share,
          se = makerobustseslist(list(fit_share, fit_cov_share)),
          title="Impact of Lawn Signs on Vote Share (Experiment 4)",
          label = "tab: exp_4",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2"),
          dep.var.labels= c("Vote Share"),
          covariate.labels=c("Assigned Lawn Signs (n=20)", "Adjacent to Lawn Signs (n=44)", NA),
          omit="(2010)|(2008)|(2006)|(2004)|(2002)|(2000)", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{10cm}{Covariates: Gubernatorial Vote Margin '02, '06, '10 and Presidential Vote Margin '00, '04, '08.}"),
          notes.append=FALSE)
#sink()

## replication table c.4 appendix
# sink("experiment_4_appendix.tex")
stargazer(fit_margin, fit_cov_margin, fit_turnout, fit_cov_turnout,
          se = makerobustseslist(list(fit_margin, fit_cov_margin, fit_turnout, fit_cov_turnout)),
          title="Impact of Lawn Signs on Margin and Turnout (Experiment 4)",
          label = "tab: exp_4_app",
          model.names=F, model.numbers=F, style="apsr",
          column.labels = c("Model 1", "Model 2", "Model 1", "Model 2"),
          dep.var.labels= c("Vote Margin", "Turnout"),
          covariate.labels=c("Assigned Lawn Signs (n=20)", "Adjacent to Lawn Signs (n=44)", NA),
          omit="(2010)|(2008)|(2006)|(2004)|(2002)|(2000)", 
          omit.labels = "Covariate Adjustment",
          omit.yes.no=c("yes", "no"),
          omit.stat=c("adj.rsq", "ser", "f"),
          star.cutoffs = c(NA, NA, NA),
          notes= c("\\parbox[t]{13cm}{Column 2 covariates: Gubernatorial Vote Share '02, '06, '10 and Presidential Vote Share '00, '04, '08. Column 4 covariates: Gubernatorial Turnout '02, '06, '10 and Turnout '00, '04, '08.}"),
          notes.append=FALSE)
# sink()


# Randomization Inference -------------------------------------------------

f_margin <-
  get_f(formula_F = "margin ~ condition_factor + GOVM2010 + USPM2008 + GOVM2006 + USPM2004 + GOVM2002 + USPM2000",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = PA_41_df)

f_share <-
  get_f(formula_F = "share ~ condition_factor + GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = PA_41_df)

f_turnout <-
  get_f(formula_F = "turnout ~ condition_factor + GOVT2010 + USPT2008 + GOVT2006 + USPT2004 + GOVT2002 + USPT2000",
        treat_var = "condition_factor", 
        weights = "weights", 
        data = PA_41_df)


sims<- 10000
f_margin_sims <- rep(NA, sims)
f_share_sims <- rep(NA, sims)
f_turnout_sims <- rep(NA, sims)

for (i in 1:sims){
  #print(i) 
  condition.sims <- restricted_condition_mat[,i]
  
  weights.sims <- 1/get_prob_obs(condition.sims)
  
  condition_factor_sims <- factor(condition.sims, levels=c("control", "adjacent", "treated"), labels = c("Control", "Adjacent", "Treated"))
  condition_factor_sims <- factor(condition_factor_sims, levels=c("Control", "Treated", "Adjacent"))
  
  frame_sim <- data.frame(PA_41_df, condition_factor_sims, weights.sims)
  
  f_margin_sims[i] <-
    get_f(formula_F = "margin ~ condition_factor_sims + GOVM2010 + USPM2008 + GOVM2006 + USPM2004 + GOVM2002 + USPM2000",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = frame_sim)
  
  f_share_sims[i] <-
    get_f(formula_F = "share ~ condition_factor_sims + GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = frame_sim)
  
  f_turnout_sims[i] <-
    get_f(formula_F = "turnout ~ condition_factor_sims + GOVT2010 + USPT2008 + GOVT2006 + USPT2004 + GOVT2002 + USPT2000",
          treat_var = "condition_factor_sims", 
          weights = "weights.sims", 
          data = frame_sim)
}
pval_share<-mean(f_share_sims >= f_share)
pval_margin <-mean(f_margin_sims >= f_margin)
pval_turnout <-mean(f_turnout_sims >= f_turnout)

ps <- t(matrix(c(pval_margin, pval_share, pval_turnout)))
colnames(ps) <- c("Margin", "Share",  "Turnout")
rownames(ps) <- "p(f_sims >= f_obs)"

# sink(file ="experiment_4_ps.tex")
## replication randomization pvalue referred to in text, section 5.4
xtable(ps, caption = "Experiment 4 Randomization Inference p-values")
# sink()