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


ncores = 12 # how much time on 12 nodes?

# note to Sayali and Bruce, all files should be in the top-level
# SPPQReplicationArchive directory---no subfolders

# Butler & Nickerson replication
# what data files are needed?
source("butler_nickerson_analysis.R")
# what files are produced?




####
# Bergan and Cole replication
####

# what data files are needed?
source("bergan_cole_analysis.R")
# what files are produced?

# Separate script to produce the tables
# what files are needed?
source("tables.R")
# table4.txt
# table5.txt

