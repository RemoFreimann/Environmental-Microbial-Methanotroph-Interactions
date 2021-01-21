

library(igraph) 
#set the size of random network
n=206   #number of nodes
e=2400  #number of edges

#generate 1000 random networks
for (i in 1:1000) {
  g <- erdos.renyi.game(n, e,'gnm',weight=T,mode="undirected")
  
  # Global toplogical features
  c <- cluster_walktrap(g)
  md <- modularity(g, membership(c), weights = NULL)
  cc <- transitivity(g, vids = NULL,
                     weights = NULL)
  spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
  gd  <- graph.density(g, loops=FALSE)
  nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
  
  ND <- degree(g, v = V(g), mode="all")
  ad  <- mean(ND)
  
  global.topol <- data.frame(n,e,cc,spl,md,gd,nd,ad)
  
  write.table(global.topol, file = sprintf("N%dE%d.er.random.network.xls",n,e), 
              append = TRUE, sep = "\t",row.names = FALSE, col.names = TRUE) }

# print node distribution statistics

degree <- data.frame(table(degree=factor(ND, levels=seq_len(max(ND)))))
degree$degree <- as.numeric(degree$degree)
plot(degree)

degree.gephi <- read.csv("edges_diablo_spearman_degree.csv")
network.degree.freq.gephi <- data.frame(table(degree=factor(degree.add, levels=seq_len(max(degree.add)))))
network.degree.freq.gephi$degree <- as.numeric(network.degree.freq.gephi$degree)
plot(network.degree.freq.gephi, log="xy", type="none", xlim=c(1,130), ylim=c(1,60) )
points(network.degree.freq.gephi, add=TRUE)
points(degree, add=TRUE, col="green")


ggplot(data=network.degree.freq.gephi, aes(x=degree, y=Freq)) +
  geom_point(col="darkgreen", cex=4) + 
  scale_x_log10(limits = c(1, 100)) + scale_y_log10(limits = c(1, 100)) +
  geom_smooth(method="auto") 


ggplot(data=degree, aes(x=degree, y=Freq)) +
  geom_point(col="darkgreen", cex=4) +
  scale_x_log10(limits = c(1, 100)) + scale_y_log10(limits = c(1, 100)) +
   geom_smooth(method="auto") 

