This is a replication archive for Phadke and Desmarais, "Considering Network Effects in the Design and Analysis of Field Experiments on State Legislatures" (2018). This subfolder contains files for replicating the Coppock analyses presented in table 4.


Instructions:
1. Run pbs files in list B on a cluster. All of them are currently set to parallelize across 12 or 20 cores, depending on the analysis.
2. Run R scripts in list D to summarize results seen in table 5 in the paper.


A. List of data files: Many of these are files that we gathered from original sources, including from other studies/authors. Some of these are created by authors of this paper. These files are needed to run the analyses in list B.
A1. nm.replication.tab
A2. CoppockJEPS.rdata
A3. committee.number.shared.csv
A4. cohort_copart_network.RData
A5. cohort_copart_network_weighted.RData


B. List of files that need to be run to produce intermediate results. These need to be run to produce the intermediate results files in list C. All of these can be run separately. One would need to specify a working directory for each .pbs file. Currently it refers to the repository used by the authors. A command of type “#PBS -M johndoe@mail.com” can be added to each .pbs file to receive updates on the code.
B1. coppock_copartisan_committee_binary.pbs
B2. coppock_copartisan_committee_weighted.pbs
B3. coppock_replication_binary_coparty_cohort.pbs
B4. coppock_replication_weighted_coparty_cohort.pbs


C. List of files produced by the analyses in list B. These files are needed to run the results interpretation scripts in list D.
C1. CoppockSPPQRRresults_copartisan_committee_binary.RData
C2. CoppockSPPQRRresults_copartisan_committee_weighted.RData
C3. CoppockSPPQRRresults_copartisan_cohort_binary.RData
C4. CoppockSPPQRRresults_copartisan_cohort_weighted.RData


D. List of files to run on the intermediate results listed in C.
D1. coppockSPPQresults_summary.R


E. List of text files containing the tables that appear in the manuscript, produced by the script files in list D.
E1. table4_topleft_coppock_committee_binary.txt
E2. table4_topright_coppock_committee_weighted.txt
E3. table4_bottomleft_coppock_cohort_binary.txt
E4. table4_bottomright_coppock_cohort_weighted.txt



