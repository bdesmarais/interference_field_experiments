set.seed(123)
y.isolated <- 10*(runif(9)>0.5)

treated <- c(rep(0,4),1,rep(0,4))

amat <- matrix(0,9,9)

amat[1,2] <- 1
amat[3,4] <- 1
amat[2,4] <- 1
amat[2,5] <- 1
amat[4,5] <- 1
amat[6,8] <- 1
amat[5,6] <- 1
amat[5,8] <- 1
amat[6,7] <- 1
amat[8,9] <- 1

amat <- amat + t(amat)

library(sna)

gplot(amat)

coords <- cbind(c(1,2,1,2,3,4,5,4,5),c(1,1,2,2,1.5,1,1,2,2))

gplot(amat,coord=coords,jitter=F,gmode="graph",vertex.cex=2,edge.lwd=3)

# quantify the degree to which the outcome is correlated with exposrure via the experiment in the way we hypothesize via the model, even after removing the hypothesized effect of the treatment.

# test statistic abs(t.test(isolated_from_treatment,treated))
# will be small if we guessed the right beta
# high if we guessed the wrong beta
# p value is the proportion of statistics under other treatments that are greater than the one under the true treatment. 

# 
treatment <- 10
exposure <- treatment/2

observed.exposed <- c(amat%*%treated)

y.observed <- y.isolated + treatment*treated+exposure*observed.exposed


#-5,-2.5,2.5,5

h.exposure <- 15
h.y.isolated <- y.observed-treatment*treated-h.exposure*observed.exposed
obs.test <- abs(t.test(h.y.isolated[which(observed.exposed==1 & treated==0)],h.y.isolated[which(observed.exposed==0 & treated==0)],var.equal = TRUE)$statistic)

# loop over different treatment regimes
sim.tests <- numeric(9)
for(i in 1:9){
	sim.treated <- rep(0,9)
	sim.treated[i] <- 1
	sim.exposed <- c(amat%*%sim.treated)
	sim.isolated <- y.observed-treatment*sim.treated-h.exposure*sim.exposed
	sim.tests[i] <- abs(t.test(sim.isolated[which(sim.exposed==1 & sim.treated==0)],sim.isolated[which(sim.exposed==0 & sim.treated==0)],var.equal = TRUE)$statistic) 
}

# p-value is the proportion of simulated tests that are lower than the observed test

p.val <- mean(sim.tests > obs.test)








