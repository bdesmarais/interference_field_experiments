#### Public Opinion Replication
#### Butler, Daniel M., and David W. Nickerson. 2011. 
#### “Can Learning Constituency Opinion Affect How Legislators Vote? Results from a Field Experiment.” 
#### Quarterly Journal of Political Science 6(1): 55–83. DOI: 10.1561/100.00011019

setwd("~/Dropbox/professional/Research/Active/CausalityInNetworks-Agenda/Interference_in_Field_Experiments/Analysis/coppock_replication_data/") # Bruce's
setwd("D:/Dropbox/Interference_in_Field_Experiments/Analysis/coppock_replication_data/") # Sayali's

rm(list=ls())

library(foreign)
library(maptools)
library(spdep)
library(wnominate)

#####################################################
# 1. Bring in Data from Butler and Nickerson (2012) #
# 2. Bring in Data from NM House Roll Calls in 2008 #
# 3. Replication Results from Butler and Nickerson  #
# 4. Create Exposure Model                          #
# 5. Randomization Inteference                      #
#####################################################

#####################################################
# 1. Bring in Data from Butler and Nickerson (2012) #
#####################################################
butler <- read.delim("nm.replication.tab")

butler <- within(butler,{
#Generate IDs
id <- 1:70
  
#Generate Presidential 2-party Vote Share
dem_vs <- kerry/(bush+kerry)

#Generate Governor Vote
v_richardson <- richardsongovdem/(richardsongovdem+dendahlgovrep)

#Generate dichotomous spending preferences.  full_co =1 if q2full is below the median. 
#(b/c the proposal was popular, the "treatment" should 
#have a larger effect among legislators 
#who learned that public support was low)

full_co <- q2full <= quantile(q2full, probs=.50)
lowsupport <- full_co
treatment_lowsupport <- treatment*lowsupport
treatment_highsupport <- treatment*(1-lowsupport)
})

#####################################################
# 3. Generate Similarity Scores                     #
#####################################################

nmhouse2008 <-read.csv("CoppockJEPS_rollcalldata.csv")
bills <- data.frame(nmhouse2008[5:21])

####Nominate Scores
bills_nona <- bills
bills_nona[bills_nona==99] <- NA
rollcalls <- rollcall(bills_nona)
nominate_scores <- wnominate(rollcalls, polarity=c(1, 2), minvotes=10)
dwnom_scores <- nominate_scores$legislators$coord1D

#pdf("CoppockJEPS_appendixfigure1.pdf")
#hist(dwnom_scores, breaks=50, main="", xlab="W-NOMINATE score")
#dev.off()

get.similarity <- function(x, y){
  return((2-abs(x-y))/2)
}

similarity.matrix.dw <- matrix(NA, ncol=70, nrow=70)
for (i in 1:70){
  for (j in 1:70){
    similarity.matrix.dw[i,j] <- get.similarity(dwnom_scores[i], dwnom_scores[j])
  }
}
diag(similarity.matrix.dw) <- 0
similarity.matrix.dw[is.na(similarity.matrix.dw)==T] <- 0


### Geographical -- connected if CDs touch

### Bring in shapefile
NM_geo <-readShapeSpatial("hd_court_ordered.shp")
### Convert to Adjacency Matrix
similarity.matrix.geo <- nb2mat(poly2nb(NM_geo), style="B")
rownames(similarity.matrix.geo) <- NULL

## Read in SNA package
library(sna)

# create an indicator for voted with district
voted_with_district <- butler$sb24==(1-butler$lowsupport)

outcome_col <- rep("white",length(voted_with_district))
outcome_col[which(voted_with_district==1)] <- "black"
outcome_col[which(is.na(voted_with_district))]  <- "gray60"



setwd("~/Dropbox/professional/Research/Active/CausalityInNetworks-Agenda/Interference_in_Field_Experiments/Analysis/") # Bruce's
setwd("D:/Dropbox/Interference_in_Field_Experiments/Analysis/") # Sayali's

#pdf("coppock_geographic_net.pdf")
net <- similarity.matrix.geo

treated <- butler$treatment

neighbor_treated <- 1*((1-treated)*(net%*%cbind(treated)) > 0)

node_sides <- round(100/(.001+treated*33+neighbor_treated*25))
set.seed(5)
gplot(net,gmode="graph",vertex.col=outcome_col,vertex.sides=node_sides,edge.col=rgb(100,100,100,100,maxColorValue=256))
# triangles treated
# grey abstained
# black voted with district
# white voted against district

dev.off()

pdf("coppock_ideological_net.pdf")
## Now with ideological network, closest 5%
threshold <- quantile(similarity.matrix.dw,.95)
net <- similarity.matrix.dw >= threshold
neighbor_treated <- 1*((1-treated)*(net%*%cbind(treated)) > 0)
node_sides <- round(100/(.001+treated*33+neighbor_treated*25))
set.seed(5)
gplot(net,gmode="graph",vertex.col=outcome_col,vertex.sides=node_sides,edge.col=rgb(100,100,100,100,maxColorValue=256))
dev.off()

### Read in Committee Data
bn_named <- read.csv("~/Dropbox/professional/Research/Active/CausalityInNetworks-Agenda/Interference_in_Field_Experiments/analysis/nm_legislators/NM_new_data/nm.replication.csv",stringsAsFactors=F) #BD
bn_named <- read.csv("D:/Dropbox/Interference_in_Field_Experiments/Analysis/nm_legislators/NM_new_data/nm.replication.csv",stringsAsFactors=F)

committees <- read.csv("~/Dropbox/professional/Research/Active/CausalityInNetworks-Agenda/Interference_in_Field_Experiments/analysis/nm_legislators/NM_new_data/newmexicostandingcommittees20072008/house_committees_2008.csv",header=F,stringsAsFactors=F,row.names=1) #BD
committees <- read.csv("D:/Dropbox/Interference_in_Field_Experiments/Analysis/nm_legislators/NM_new_data/newmexicostandingcommittees20072008/house_committees_2008.csv",header=F,stringsAsFactors=F,row.names=1) #SP

committee_matches <- NULL
library(stringdist)
legislator_names <- paste(bn_named$last_name,bn_named$first_name,sep="")
legislator_names <- gsub(" ","",legislator_names)
legislator_names <- tolower(legislator_names)
committee_amat <- matrix(0,length(legislator_names),length(legislator_names))
for(i in 1:nrow(committees)){
	committeei <- c(t(committees[i,]))
	committeei <- committeei[which(nchar(committeei) > 4)]
	committeei <- tolower(committeei)
	committeei <- gsub(" ","",committeei)
	committeei <- gsub("\"","",committeei)
	committeei <- gsub("[[:punct:]]", "", committeei)
	committeei <- iconv(committeei,to="UTF-8")
	name_dist_mat <- stringdistmatrix(committeei,legislator_names)
	matches <- apply(name_dist_mat,1,which.min)
	mins <- apply(name_dist_mat,1,min)
	for(j in 2:length(matches)){
		for(h in 1:(j-1)){
			committee_amat[matches[j],matches[h]] <- committee_amat[matches[j],matches[h]] + 1
			committee_amat[matches[h],matches[j]] <- committee_amat[matches[h],matches[j]] + 1
		}
	}
	committee_matches <- rbind(committee_matches,cbind(committeei,legislator_names[matches],mins))
}

#pdf("nm_committee_net.pdf")
## Now with committee network, tied if on more than 1 together
net <- committee_amat > 1
neighbor_treated <- 1*((1-treated)*(net%*%cbind(treated)) > 0)
node_sides <- round(100/(.001+treated*33+neighbor_treated*25))
set.seed(5)
gplot(net,gmode="graph",vertex.col=outcome_col,vertex.sides=node_sides,edge.col=rgb(100,100,100,100,maxColorValue=256))
dev.off()








