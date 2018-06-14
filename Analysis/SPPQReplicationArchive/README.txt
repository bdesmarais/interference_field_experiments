This is a replication archive for Phadke and Desmarais, "Considering Network Effects in the Design and Analysis of Field
Experiments on State Legislatures" (2018). There are a total of eight analyses.


Instructions:
1. Run pbs files in list B on a cluster. All of them are currently set to parallelize across 20 cores.
2. Run R scripts in list D to summarize results seen in tables 4 and 5 in the paper.


A. List of data files: Many of these are files that we gathered from original sources, including from other
studies/authors. Some of these are created by authors of this paper. These files are needed to run the analyses in list B.
A1. bergan_data.csv
A2. bergan.dta
A3. cohort_copart_network.RData
A4. cohort_network.RData
A5. cosponsorship_network.RData
A6. w_cohort_copart_network.RData
A7. w_cohort_network.RData


B. List of files that need to be run to produce intermediate results. These need to be run to produce the intermediate
results files in list C. All of these can be run separately.
B1. bergan_cohort_bin_newmodel_w_chamber1.pbs
B2. bergan_cohort_bin_newmodel_w_chamber2.pbs
B3. bergan_cohort_bin_newmodel_w_chamber3.pbs
B4. bergan_cohort_bin_newmodel_w_chamber4.pbs
B5. bergan_cohort_bin_newmodel_w_chamber5.pbs
B6. bergan_cohort_bin_newmodel_w_chamber6.pbs
B7. bergan_cohort_bin_newmodel_w_chamber7.pbs
B8. bergan_cohort_bin_newmodel_w_chamber8.pbs
B9. bergan_cohort_bin_newmodel_w_chamber9.pbs
B10. bergan_cohort_bin_newmodel_w_chamber10.pbs
B11. bergan_cohort_bin_newmodel_w_chamber11.pbs
B12. bergan_cohort_bin_newmodel_w_chamber12.pbs
B13. bergan_cohort_bin_newmodel_w_chamber13.pbs
B14. bergan_cohort_bin_newmodel_w_chamber14.pbs
B15. bergan_cohort_bin_newmodel_w_chamber15.pbs
B16. bergan_cohort_bin_newmodel_w_chamber16.pbs
B17. bergan_cohort_sim_newmodel_w_chamber1.pbs
B18. bergan_cohort_sim_newmodel_w_chamber2.pbs
B19. bergan_cohort_sim_newmodel_w_chamber3.pbs
B20. bergan_cohort_sim_newmodel_w_chamber4.pbs
B21. bergan_cohort_sim_newmodel_w_chamber5.pbs
B22. bergan_cohort_sim_newmodel_w_chamber6.pbs
B23. bergan_cohort_sim_newmodel_w_chamber7.pbs
B24. bergan_cohort_sim_newmodel_w_chamber8.pbs
B25. bergan_cohort_sim_newmodel_w_chamber9.pbs
B26. bergan_cohort_sim_newmodel_w_chamber10.pbs
B27. bergan_cohort_sim_newmodel_w_chamber11.pbs
B28. bergan_cohort_sim_newmodel_w_chamber12.pbs
B29. bergan_cohort_sim_newmodel_w_chamber13.pbs
B30. bergan_cohort_sim_newmodel_w_chamber14.pbs
B31. bergan_cohort_sim_newmodel_w_chamber15.pbs
B32. bergan_cohort_sim_newmodel_w_chamber16.pbs
B33. bergan_cospon_binary_newmodel1.pbs
B34. bergan_cospon_binary_newmodel2.pbs
B35. bergan_cospon_binary_newmodel3.pbs
B36. bergan_cospon_binary_newmodel4.pbs
B37. bergan_cospon_binary_newmodel5.pbs
B38. bergan_cospon_binary_newmodel6.pbs
B39. bergan_cospon_binary_newmodel7.pbs
B40. bergan_cospon_binary_newmodel8.pbs
B41. bergan_cospon_binary_newmodel9.pbs
B42. bergan_cospon_binary_newmodel10.pbs
B43. bergan_cospon_binary_newmodel11.pbs
B44. bergan_cospon_binary_newmodel12.pbs
B45. bergan_cospon_binary_newmodel13.pbs
B46. bergan_cospon_binary_newmodel14.pbs
B47. bergan_cospon_binary_newmodel15.pbs
B48. bergan_cospon_binary_newmodel16.pbs
B49. bergan_cospon_weighted_newmodel.pbs


C. List of files produced by the analyses in list B. These files are needed to run the results interpretation scripts in
list D.
C1. BerganSPPQRRresults_copartisan_cohort_chamber_binary1.RData
C2. BerganSPPQRRresults_copartisan_cohort_chamber_binary2.RData
C3. BerganSPPQRRresults_copartisan_cohort_chamber_binary3.RData
C4. BerganSPPQRRresults_copartisan_cohort_chamber_binary4.RData
C5. BerganSPPQRRresults_copartisan_cohort_chamber_binary5.RData
C6. BerganSPPQRRresults_copartisan_cohort_chamber_binary6.RData
C7. BerganSPPQRRresults_copartisan_cohort_chamber_binary7.RData
C8. BerganSPPQRRresults_copartisan_cohort_chamber_binary8.RData
C9. BerganSPPQRRresults_copartisan_cohort_chamber_binary9.RData
C10. BerganSPPQRRresults_copartisan_cohort_chamber_binary10.RData
C11. BerganSPPQRRresults_copartisan_cohort_chamber_binary11.RData
C12. BerganSPPQRRresults_copartisan_cohort_chamber_binary12.RData
C13. BerganSPPQRRresults_copartisan_cohort_chamber_binary13.RData
C14. BerganSPPQRRresults_copartisan_cohort_chamber_binary14.RData
C15. BerganSPPQRRresults_copartisan_cohort_chamber_binary15.RData
C16. BerganSPPQRRresults_copartisan_cohort_chamber_binary16.RData
C17. BerganSPPQRRresults_copartisan_cohort_chamber_similarity1.RData
C18. BerganSPPQRRresults_copartisan_cohort_chamber_similarity2.RData
C19. BerganSPPQRRresults_copartisan_cohort_chamber_similarity3.RData
C20. BerganSPPQRRresults_copartisan_cohort_chamber_similarity4.RData
C21. BerganSPPQRRresults_copartisan_cohort_chamber_similarity5.RData
C22. BerganSPPQRRresults_copartisan_cohort_chamber_similarity6.RData
C23. BerganSPPQRRresults_copartisan_cohort_chamber_similarity7.RData
C24. BerganSPPQRRresults_copartisan_cohort_chamber_similarity8.RData
C25. BerganSPPQRRresults_copartisan_cohort_chamber_similarity9.RData
C26. BerganSPPQRRresults_copartisan_cohort_chamber_similarity10.RData
C27. BerganSPPQRRresults_copartisan_cohort_chamber_similarity11.RData
C28. BerganSPPQRRresults_copartisan_cohort_chamber_similarity12.RData
C29. BerganSPPQRRresults_copartisan_cohort_chamber_similarity13.RData
C30. BerganSPPQRRresults_copartisan_cohort_chamber_similarity14.RData
C31. BerganSPPQRRresults_copartisan_cohort_chamber_similarity15.RData
C32. BerganSPPQRRresults_copartisan_cohort_chamber_similarity16.RData
C33. BerganSPPQRRresults_cospon_binary1.RData
C34. BerganSPPQRRresults_cospon_binary2.RData
C35. BerganSPPQRRresults_cospon_binary3.RData
C36. BerganSPPQRRresults_cospon_binary4.RData
C37. BerganSPPQRRresults_cospon_binary5.RData
C38. BerganSPPQRRresults_cospon_binary6.RData
C39. BerganSPPQRRresults_cospon_binary7.RData
C40. BerganSPPQRRresults_cospon_binary8.RData
C41. BerganSPPQRRresults_cospon_binary9.RData
C42. BerganSPPQRRresults_cospon_binary10.RData
C43. BerganSPPQRRresults_cospon_binary11.RData
C44. BerganSPPQRRresults_cospon_binary12.RData
C45. BerganSPPQRRresults_cospon_binary13.RData
C46. BerganSPPQRRresults_cospon_binary14.RData
C47. BerganSPPQRRresults_cospon_binary15.RData
C48. BerganSPPQRRresults_cospon_binary16.RData
C49. BerganSPPQRRresults_cospon_weighted_newmodel.RData


D. List of files to run on the intermediate results listed in C.
D1. berganSPPQresults_summary_cohortBinary.R
D2. berganSPPQresults_summary_cohortWeighted.R
D3. berganSPPQresults_summary_cosponBinary.R
D4. berganSPPQresults_summary.R


E. List of text files containing the tables that appear in the manuscript, produced by the script files in list D.