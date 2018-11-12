##################
### Setting Up ###
##################

library(igraph)
library(MASS)

color.groups <- function(z, col = tim.colors(50),
zlim = c(min(z, na.rm = TRUE),
max(z, na.rm = TRUE))){
    breaks <- seq(zlim[1], zlim[2], length = length(col) + 1)
    zgroups <- cut(z, breaks = breaks, labels = FALSE,
    include.lowest = TRUE)
    return(col[zgroups])
}

n = 20

##############
### Stop 1 ###
##############

### Erdos-Renyi
### Let's take an Erdos Renyi random graph as a running example

set.seed(1000)

p = 0.12
adjmat.erdosrenyi = matrix(0,n,n)
adjmat.erdosrenyi[upper.tri(adjmat.erdosrenyi)] = rbinom(choose(n,2),1,p)
adjmat.erdosrenyi = adjmat.erdosrenyi + t(adjmat.erdosrenyi)

g.erdosrenyi = graph_from_adjacency_matrix(adjmat.erdosrenyi, mode = "undirected")

plot(g.erdosrenyi,layout=layout.kamada.kawai,vertex.color=blues9[5])

##############
### Stop 2 ###
##############

### Eigenvector Centrality

eigen.erdosrenyi = eigen(adjmat.erdosrenyi)

niters = 30
x = matrix(0,n,niters)
x[,1] = rep(1,n)
for(i in 2:niters){x[,i] = 1/eigen.erdosrenyi$values[1] * adjmat.erdosrenyi %*% x[,i-1]}

co <- layout_with_kk(g.erdosrenyi)

plotpower = function(ith){
for(ith in 1:30){
	plot(g.erdosrenyi, layout=co, vertex.size=30, rescale=FALSE,
		xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
		vertex.color = ifelse(rep(ith==1,n),blues9[4],color.groups(x[,ith],blues9)))
	title(main=paste("x",ith-1,sep=""),line=-0.1)
	readline()
}
}

##############
### Stop 3 ###
##############

plotpower() # <------ wait here, sequential plots

##############
### Stop 4 ###
##############

v = sum(eigen.erdosrenyi$vectors[,1]*x[,1])*eigen.erdosrenyi$vectors[,1]
plot(g.erdosrenyi, layout=co, vertex.size=30, rescale=FALSE,
     xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
     vertex.color = color.groups(v,blues9))
title(main="x_inf, theoretically",line=-0.1)

##############
### Stop 5 ###
##############

### Random DAG

set.seed(1000)

p = 0.12
adjmat.dag = matrix(0,n+5,n+5)
adjmat.dag[upper.tri(adjmat.dag)] = c(rbinom(choose(n,2),1,p),rep(0,choose(n+5,2)-choose(n,2)))
adjmat.dag[matrix(c(rep(n,5),n+1:5),5,2)] = 1
adjmat.dag[20,1] = 1

g.dag = graph_from_adjacency_matrix(adjmat.dag, mode = "directed")

plot(g.dag,layout=layout.kamada.kawai,vertex.color=blues9[5],edge.arrow.size=0.5)

##############
### Stop 6 ###
##############

### Eigenvector Centrality

eigen.dag = eigen(adjmat.dag)

niters = 30
x = matrix(0,n+5,niters)
x[,1] = rep(1,n+5)
for(i in 2:niters){x[,i] = t(adjmat.dag) %*% x[,i-1]}

co <- layout_with_kk(g.dag)

plotpowerpadded = function(){
for(ith in 1:30){
	plot(g.dag, layout=co, vertex.size=30, rescale=FALSE,
		xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
		vertex.color = ifelse(rep(ith==1,n+5),blues9[4],
				      ifelse(rep(all(x[,ith]==0),n+5),"#FFFFFF",
					     color.groups(x[,ith],blues9))),
	     edge.arrow.size=0.5)
	title(main=paste("x",ith-1,sep=""),line=-0.1)
	readline()
}
}

##############
### Stop 7 ###
##############

plotpowerpadded() # <------ sequential plots

##############
### Stop 8 ###
##############

### Katz Centrality

niters = 30
x = matrix(0,n+5,niters)
x[,1] = rep(1,n+5)
for(i in 2:niters){x[,i] = 0.3 * t(adjmat.dag) %*% x[,i-1] + rep(0.1,n+5)}

co <- layout_with_kk(g.dag)

##############
### Stop 9 ###
##############

plotpowerpadded() # <------ sequential plots

###############
### Stop 10 ###
###############

v = ginv(diag(n+5) - 0.3 * t(adjmat.dag)) %*% rep(0.1,n+5)
plot(g.dag, layout=co, vertex.size=30, rescale=FALSE,
	xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
	vertex.color = color.groups(v,blues9),
     edge.arrow.size=0.5)
title(main="x_inf",line=-0.1)

###############
### Stop 11 ###
###############

### Page Rank

niters = 30
x = matrix(0,n+5,niters)
x[,1] = rep(1,n+5)
DOut =  apply(adjmat.dag,1,sum)
DOut[DOut==0] = 1
for(i in 2:niters){x[,i] = 0.3 * t(adjmat.dag) %*% diag(1/DOut) %*% x[,i-1] + rep(0.1,n+5)}

co <- layout_with_kk(g.dag)

###############
### Stop 12 ###
###############

plotpowerpadded() # <------ sequential plots

###############
### Stop 13 ###
###############

v = ginv(diag(n+5) - 0.3 * t(adjmat.dag) %*% diag(1/DOut)) %*% rep(0.1,n+5)
plot(g.dag, layout=co, vertex.size=30, rescale=FALSE,
	xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
	vertex.color = color.groups(v,blues9),
     edge.arrow.size=0.5)
title(main="x_inf",line=-0.1)


