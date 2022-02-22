library('igraph')
dir ('Output/',pattern ='csv')
#[1] "Psi.csv"     "Psi_b.csv"   "Psi_n.csv"   "Theta.csv"   "Theta_b.csv" "Theta_n.csv"
filename='./Output/Psi_n.csv'
M=read.csv(filename,header=FALSE)
M=as.matrix(M)
G=graph_from_adjacency_matrix(M,mode = "undirected")
V=is.connected(G)


lp=cluster_label_prop(G)
ComSizes=sizes(lp)
ComSizes[which(sizes(lp)>1)]
mc=sort(c(which(lp$membership==69), which(lp$membership==70)))
#mc=sort(c(which(lp$membership==69), which(lp$membership==70), which(lp$membership==78)))
length(intersect(c(118:182),mc))/65#length(mc)
       


modularity(lp)

edge_density(G)

pdf('./plots_luisa/Psi_nlp.pdf')
layout <-layout.fruchterman.reingold(G)
plot(lp, G, layout=layout, vertex.label=NA, vertex.size=5,  edge.arrow.size=.2)
dev.off()






filename='./Output/Theta_n.csv'
M=read.csv(filename,header=FALSE)
M=as.matrix(M)
G=graph_from_adjacency_matrix(M,mode = "undirected")
V=is.connected(G)
LOU=cluster_louvain(G)
max(LOU$membership)
FG=cluster_fast_greedy(G)
max(FG$membership)
ComSizes=sizes(LOU)
ComSizes[which(sizes(LOU)>1)]



lp=cluster_label_prop(G)
ComSizes=sizes(lp)
ComSizes[which(sizes(lp)>1)]

modularity(lp)

edge_density(G)


pdf('./plots_luisa/Theta_nlp.pdf')
layout <-layout.fruchterman.reingold(G)
plot(lp, G, layout=layout, vertex.label=NA, vertex.size=5,  edge.arrow.size=.2)
dev.off()

