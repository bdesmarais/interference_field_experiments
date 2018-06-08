## This file calls source code to replicate the figures and tables
## in Phadke and Desmarais, "Considering Network Effects in the Design and Analysis of 
## Field Experiments on State Legislatures"
## in SPPQ. For each call to source(), we list the required input files above in comments
## and the output files produced by the code below. Note, this code takes a long time to run,
## and takes an argument (ncores), which can be used to parallelize across cores on one node.

## Required packages
#Packages
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages("iterators", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("foreign", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("foreach", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("doParallel", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/")
install.packages("xtable", Sys.getenv("R_LIBS_USER"), repos = "https://cran.cnr.berkeley.edu/") #Producing tables for latex
library(foreign,lib.loc=Sys.getenv("R_LIBS_USER"))
library(foreach,lib.loc=Sys.getenv("R_LIBS_USER"))
library(doParallel,lib.loc=Sys.getenv("R_LIBS_USER"))
library(xtable,lib.loc=Sys.getenv("R_LIBS_USER"))


# #### Cleanup
# rm(list=ls())
# gc()
# set.seed(231)
# 
# ncores = 40 # how much time on 40 nodes?
# 
# # note to Sayali and Bruce, all files should be in the top-level
# # SPPQReplicationArchive directory---no subfolders
# 
# ####
# # Butler & Nickerson replication (Coppock)
# ####
# 
# # what data files are needed?
# # 'nm.replication.tab', 'CoppockJEPS.RData', 'butler_cohort_copart_network.RData', 'butler_cohort_copart_network_weighted.RData', 'butler_committee_number_shared.csv'
# source("butler_nickerson_analysis.R")
# # what files are produced?
# # 'CoppockSPPQRRresults_copartisan_cohort_binary.RData', 'CoppockSPPQRRresults_copartisan_cohort_weighted.RData', 'CoppockSPPQRRresults_copartisan_committee_binary.RData', 'CoppockSPPQRRresults_copartisan_committee_weighted.RData'
# 
# # Separate script to produce the tables
# # what files are needed?
# # the four mentioned above
# source("butler_tables.R")
# # Outputs four tables, one for each analysis
# # 'butler_cohort_binary_table4_bottomleft.txt', 'butler_cohort_weighted_table4_bottomright.txt', 'butler_committee_binary_table4_topleft.txt', 'butler_committee_weighted_table4_topright.txt'


#### Cleanup
rm(list=ls())
gc()
set.seed(132)

ncores = 40 # how much time on 40 nodes?

####
# Bergan and Cole replication
####

# what data files are needed?
# 'bergan.dta', 'bergan_cohort_copart_network.RData', 'bergan_w_cohort_copart_network.RData', 'bergan_cosponsorship_network.RData'
source("bergan_cole_analysis.R")
# what files are produced?
# 'BerganSPPQRRresults_copartisan_cohort_binary.RData', 'BerganSPPQRRresults_copartisan_cohort_similarity.RData', 'BerganSPPQRRresults_cospon_binary.RData', 'BerganSPPQRRresults_cospon_weighted.RData'

# Separate script to produce the tables
# what files are needed?
# the four mentioned above
source("bergan_tables.R")
# Outputs four tables, one for each analysis
# 'bergan_cohort_binary_table5_bottomleft.txt', 'bergan_cohort_weighted_table5_bottomright.txt', 'bergan_cospon_binary_table5_topleft.txt', 'bergan_cospon_weighted_table5_topright.txt' 

