# The Effects of Lawn Signs on Vote Outcomes: 
# Results from Four Randomized Field Experiments
# Donald P. Green, Jonathan S. Krasno, Alexander Coppock, Benjamin D. Farrer, Brandon Lenoir, and Joshua N. Zingher

# Experiment 4 Random Assignment
# this file creates 10,000 Resticted Randomizations and outputs them, along with probabilities of assignment.

# Uncomment to set your working directory
# setwd("")

# Uncomment to install required packages
# install.packages(c("randomizr", "dplyr", "beepr"))

rm(list=ls())

library(dplyr)
library(beepr)
library(randomizr)
# source("programs/cumberland source.R")
source("Lawn_Signs_Source.R")

load("Experiment_4_Past_Elections.rdata")
load("Experiment_4_Adjacencies.rdata")
load("Experiment_4_Exclusions.rdata")

PA_41_df <- 
  left_join(PA_41_df, PA_41_exclusions)

PA_41_df_subset <- filter(PA_41_df, cantuse==0)
N_subset <- nrow(PA_41_df_subset)

PA_41_adjmat_subset <- PA_41_adjmat[PA_41_df$cantuse==0,PA_41_df$cantuse==0]

Z_mat <- replicate(100000,
                   complete_ra(N = N_subset, m = 20))

cond_mat <- apply(Z_mat, 2, get_condition, adj=PA_41_adjmat_subset)

prob_treat <- rowMeans(cond_mat=="treated", na.rm=TRUE)
prob_adj <- rowMeans(cond_mat=="adjacent", na.rm=TRUE)
prob_control <- rowMeans(cond_mat=="control", na.rm=TRUE)

### Choose 10,000 of these 100,000 that fits the t-ratio<1 criterion.

restricted_sims <- 10000

restricted_Z_mat <- matrix(NA, nrow=N_subset, ncol=restricted_sims)
restricted_condition_mat <- matrix(NA, nrow=N_subset, ncol=restricted_sims)

i<-1
j<-1
available_Z <- rep(1, dim(Z_mat)[2])
available_Z_index <- 1:dim(Z_mat)[2]
while(i<=restricted_sims){
  selected_Z <- which(available_Z_index %in% sample(available_Z_index[available_Z==1], 1))
  available_Z[selected_Z] <-0
  
  condition_obs <- cond_mat[ , selected_Z]
  weights <- rep(NA, N_subset)
  weights[condition_obs=="treated"] <- 1/prob_treat[condition_obs=="treated"]
  weights[condition_obs=="adjacent"] <- 1/prob_adj[condition_obs=="adjacent"]
  weights[condition_obs=="control"] <- 1/prob_control[condition_obs=="control"]
  
  PA_sim <- subset(data.frame(PA_41_df_subset, condition_obs, weights, prob_treat, prob_adj, prob_control),
                   prob_treat<1 & prob_treat>0 & prob_adj<1 & prob_adj>0 & prob_control<1 & prob_control>0) 
  
  fit_3_1 <- lm(condition_obs=="treated" ~ GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000, 
                weights=weights, data=subset(PA_sim, condition_obs != "control"))
  fit_2_1 <- lm(condition_obs=="adjacent" ~ GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000, 
                weights=weights,  data=subset(PA_sim, condition_obs != "control"))
  fit_3_2 <- lm(condition_obs=="treated" ~ GOVP2010 + USPP2008 + GOVP2006 + USPP2004 + GOVP2002 + USPP2000, 
                weights=weights,  data=subset(PA_sim, condition_obs != "adjacent"))
  
  keep <- (lmp(fit_3_1) > .25 & lmp(fit_2_1) > .25 & lmp(fit_3_2) > .25)
  
  if(keep){
    restricted_Z_mat[,i] <- Z_mat[,selected_Z]
    restricted_condition_mat[,i] <- cond_mat[,selected_Z]
    i <- i+1
    print(i)
  }
  j<-j+1
}


dim(restricted_Z_mat)
dim(restricted_condition_mat)


# save(restricted_Z_mat, restricted_condition_mat, 
#      file = "Experiment_4_Permutation_Matrix.rdata")


# Conduct Actual Random Assignment ----------------------------------------

load("Experiment_4_Permutation_Matrix.rdata")

set.seed(1234567)
chosen_randomization <- sample(10000, 1)

Z_obs <- restricted_Z_mat[,chosen_randomization]
condition_obs <- restricted_condition_mat[,chosen_randomization]
res_prob_treated <- rowMeans(restricted_condition_mat=="treated")
res_prob_adj <- rowMeans(restricted_condition_mat=="adjacent")
res_prob_control <- rowMeans(restricted_condition_mat=="control")

PA_41_df_subset <- data.frame(PA_41_df_subset, Z_obs, condition_obs, res_prob_treated, res_prob_adj, res_prob_control) %>%
  select(NAME10, Z_obs, condition_obs, res_prob_treated, res_prob_adj, res_prob_control)

PA_41_df <- left_join(PA_41_df, PA_41_df_subset)


# save(PA_41_df, file= "Experiment_4_Data_with_Random_Assignment.rdata")


