# assumes working directly is initially set to project folder
setwd("./sppq_submission/images")

amat <- matrix(0,nrow=5,ncol=5)
amat[1,2] <- 1
amat[2,3] <- 1
amat[2,4] <- 1
amat[3,4] <- 1
amat[4,5] <- 1
amat[3,5] <- 1

# direct neighbors, two-hop neighbors
# effect degrade with distance
# mod1---two-hop, not degraded
nsides <- c(4,4,1000,3,4)

mod1 <- c("","","","","")

# light gray indicates modest effect
col1 <- c("black","black","white","red","black")

library(sna)

coord.nodes <- cbind(seq(1,3,length=5),c(1,1,1,1,1))

pdf(file="effect_constant_threehops.pdf",height=10,width=20)
par(mar=c(0,0,0,0))
gplot(amat,gmode="graph",label=mod1,coord=coord.nodes,label.pos=5,jitter=F,vertex.cex=3,vertex.col=col1,label.cex=2,vertex.sides=nsides,edge.lwd=2.75)
dev.off()

# light gray indicates modest effect
col2 <- c("grey85","grey50","white","red","black")

library(sna)

coord.nodes <- cbind(seq(1,3,length=5),c(1,1,1,1,1))

pdf(file="effect_decay_threehops.pdf",height=10,width=20)
par(mar=c(0,0,0,0))
gplot(amat,gmode="graph",label=mod1,coord=coord.nodes,label.pos=5,jitter=F,vertex.cex=3,vertex.col=col2,label.cex=2,vertex.sides=nsides,edge.lwd=2.75)
dev.off()

# light gray indicates modest effect
col3 <- c("white","black","white","red","black")

library(sna)

coord.nodes <- cbind(seq(1,3,length=5),c(1,1,1,1,1))

pdf(file="effect_constant_twohops.pdf",height=10,width=20)
par(mar=c(0,0,0,0))
gplot(amat,gmode="graph",label=mod1,coord=coord.nodes,label.pos=5,jitter=F,vertex.cex=3,vertex.col=col3,label.cex=2,vertex.sides=nsides,edge.lwd=2.75)
dev.off()

# light gray indicates modest effect
col4 <- c("white","grey50","white","red","black")

library(sna)

coord.nodes <- cbind(seq(1,3,length=5),c(1,1,1,1,1))

pdf(file="effect_decay_twohops.pdf",height=10,width=20)
par(mar=c(0,0,0,0))
gplot(amat,gmode="graph",label=mod1,coord=coord.nodes,label.pos=5,jitter=F,vertex.cex=3,vertex.col=col4,label.cex=2,vertex.sides=nsides,edge.lwd=2.75)
dev.off()




