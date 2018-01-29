library(foreign)
nm_dat <- read.dta("~/Dropbox/professional/Research/Active/causalityinnetworks-agenda/interference_in_field_experiments/Analysis/nm_legislators/qjps_11019_supp/nm.replication.dta")
write.csv(nm_dat,"~/Dropbox/professional/Research/Active/causalityinnetworks-agenda/interference_in_field_experiments/Analysis/nm_legislators/NM_new_data/nm.replication.csv",row.names=F)