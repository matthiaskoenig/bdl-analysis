##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-21
#
##############################################################################################

#---------------------------------------------
# Hierarchical clustering
#---------------------------------------------

# Perform hierarchical clustering on matrix
# The hclust function in R uses the complete linkage method for hierarchical clustering by default.
# This particular clustering method defines the cluster distance between two clusters to be the 
# maximum distance between their individual components.

# apply hirarchical clustering 
method = "ys2"
if (identical(method, "ys1")){
  cor.cluster <- cor.ys1  
}else if (identical(method, "ys2")){
  cor.cluster <- cor.ys2
}else if (identical(method, "ys3")){
  cor.cluster <- cor.ys3
}else if (identical(method, "yr1")){
  cor.cluster <- cor.yr1
}else if (identical(method, "yr2")){
  cor.cluster <- cor.yr2
}else if (identical(method, "yr3")){
  cor.cluster <- cor.yr3
}

hc <- hclust(dist(cor.cluster)) 
corrplot(cor.cluster[hc$order, hc$order], order="original", method="square", type="full",tl.cex=0.3, tl.col="black")

# plot(hc)               # plot the dendrogram 
plot(hc, hang=-1)               # plot the dendrogram 

# cut tree into clusters
Ngroups = 5
rect.hclust(hc, k=Ngroups)
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)
groups.hc.order <- groups[hc$order]

dev.off()


# Plot individual time courses in cluster
for (k in 1:Ngroups){
  fname <- sprintf("%s_cluster_%s.png", method, k)
  png(filename=sprintf("../results/cluster/%s", fname), width=3000, height=3000, res=200)
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  par(mfrow=c(N,N))
  for (name in names(g)){
    plot_single(name_A=name) 
  }
  par(mfrow=c(1,1))  
  dev.off()
}

# Plot the mean curves
# install.packages('matrixStats')
library('matrixStats')

# mean plots for clusters 
# TODO: variance to the cluster
f_normalize_centering <- function(a){
  a.norm <- (a - mean(a))/(max(a, na.rm=TRUE) - min(a, na.rm=TRUE))
  return(a.norm)
}

fname <- sprintf("%s_cluster_overview.png", method)
png(filename=sprintf("../results/cluster/%s", fname), width=1600, height=1600, res=200)
par(mfrow=c(ceiling(sqrt(Ngroups)),ceiling(sqrt(Ngroups))))
# steps = 1:8
for (k in 1:Ngroups){
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  dgroup <- BDLmean[names(g)]

  # centralize and normalize columns, i.e. the individual factors for comparison
  dgroup.norm <- apply(dgroup, 2, f_normalize_centering)
  
  # mean and sd for timepoints 
  g.mean <- rowMeans(dgroup.norm)
  g.sd <- rowSds(dgroup.norm)   # apply(dgroup.norm, 2, sd)
  
  # plot sd range
  plot(1, type="n", xlab="", ylab="", xlim=c(1, 8), ylim=c(-1, 1), main=sprintf("%s : Cluster %s", method, k))
  polygon(c(steps, rev(steps)), c(g.mean+g.sd, rev(g.mean-g.sd)),
          col = rgb(0.5,0.5,0.5,0.5), border = NA)
  
  # individual data
  for (name in names(g)){
    points(steps, dgroup.norm[, name], pch=16, col="black")
    lines(steps, dgroup.norm[, name], col=rgb(0.5,0.5,0.5, 0.8), lwd=1)
  }
  # mean over factors in cluster
  lines(steps, g.mean, col="blue", lwd=2)
}
par(mfrow=c(1,1))
dev.off()



# -----------------------------------------------------------------------------
# Dendrograms

# install.packages("dendextend")
library(dendextend)

# better dendrogram
dend1 <- as.dendrogram(hc)
# Get the package:
cutree(dend1,h=70) # it now works on a dendrogram
# It is like using:
dendextend:::cutree.dendrogram(dend1,h=70)
# dend1 <- color_branches(dend1, k = 7)
dend1 <- color_labels(dend1, k=Ngroups)
plot(dend1)



# TODO: create the list of all pairwise correlations
# TODO: GO annotations of clusters


f_cor_pair_plot("Cebpa", "Cyp24a1")
  