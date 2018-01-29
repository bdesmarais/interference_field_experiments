####
## Creating dummy network for poster
####

rm(list=ls())
library(igraph)

# Creating adjacency matrix
adj.mat <- matrix(0, 10, 10)
adj.mat[1,2] <- adj.mat[2,1] <- 1

adj.mat[2,3] <- adj.mat[3,2] <- 1
adj.mat[2,4] <- adj.mat[4,2] <- 1
adj.mat[3,4] <- adj.mat[4,3] <- 1
adj.mat[4,5] <- adj.mat[5,4] <- 1
adj.mat[4,6] <- adj.mat[6,4] <- 1
adj.mat[5,6] <- adj.mat[6,5] <- 1
adj.mat[5,8] <- adj.mat[8,5] <- 1
adj.mat[5,9] <- adj.mat[9,5] <- 1
adj.mat[6,7] <- adj.mat[7,6] <- 1
adj.mat[7,8] <- adj.mat[8,7] <- 1
adj.mat[9,10] <- adj.mat[10,9] <- 1

# Assigning treatment/control to nodes as an attribute
nodes <- as.data.frame(c("C","C","T","C","T","T","T","C","T","C"))

# Defining some reference vectors for plotting
vertcols <- c("dodgerblue4","white")
labcols <- c("white","dodgerblue4")
categories <- c("T", "C")

# Plotting
net <- graph.adjacency(adj.mat, mode = "undirected")

cols <- vertcols[match(as.character(nodes[,1]),categories)]
lcols <- labcols[match(as.character(nodes[,1]),categories)]

set.seed(5)
plot(net, vertex.color = cols, vertex.size = 15,
     vertex.label = as.character(nodes[,1]),
     vertex.label.color = lcols, vertex.label.cex = .85,
     edge.color = "black",
     layout=layout.fruchterman.reingold)
box(which = "plot", col = "dodgerblue4")

# Below commands not used because border wasn't saving in the output
# Manually saved a pdf named Dummy_network.pdf
pdf("Dummy_network_plot.pdf")
plot(net,vertex.color = cols, vertex.size = 15, vertex.label = as.character(nodes[,1]),
     vertex.label.color = lcols, vertex.label.cex = .85,
     edge.color = "black",
     layout=layout.fruchterman.reingold)
box(which = "figure", col = "black")
dev.off()

