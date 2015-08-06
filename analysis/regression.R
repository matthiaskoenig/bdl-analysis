##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-21
#
#

#
##############################################################################################

# --- YR1, YS1, YR2, YS2 correlation ----------
# Calculation of time-course based correlation measurements, namely ys1, ys2, yr1, yr2
# and the respective adaptions to the underlying datasets, i.e. using multiple
# repeat data with every time point measurment coming from an individual
# sample.
# In ys1, ys2, yr1, yr2 all calculations are performed on the mean time course data.
# The classical correlation components are replaced with the correlations calculated
# on the individual data points.

# calculate ys1, yr1 on mean data, i.e. correlation part (S*), slope part (A) and min/max part (M) are
# all calculated on the mean data of all repeats.
w <- list(w1=0.5, w2=0.25, w3=0.25)
ys1.mean <- ys1.df(BDLmean.fil, BDLmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1.mean <- 2*(ys1.mean$value - 0.5)  # scaling to interval [-1, 1]

ys2.mean <- ys2.df(BDLmean, BDLmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2.mean <- 2*(ys2.mean$value - 0.5)  # scaling to interval [-1, 1]  



# Pearson & spearman correlation on full dataset as replacement for the 
# mean Pearson/Spearman in ys1, ys2, yr1, yr2
# Now calculate the scores with full correlation
cor.S_star <- ( cor(data, method="spearman", use="pairwise.complete.obs") + 1 )/2
cor.R_star <- ( cor(data, method="pearson", use="pairwise.complete.obs") + 1 )/2
cor.S_star.fil <- ( cor(data.fil, method="spearman", use="pairwise.complete.obs") + 1 )/2
cor.R_star.fil <- ( cor(data.fil, method="pearson", use="pairwise.complete.obs") + 1 )/2


w <- list(w1=0.5, w2=0.25, w3=0.25)
cor.ys1_full <- w$w1*cor.S_star + w$w2*res.ys1$A + w$w3*res.ys1$M
cor.ys1_full.scaled <- 2*(cor.ys1_full-0.5)
cor.ys2_full <- w$w1*cor.S_star + w$w2*res.ys2$A_star + w$w3*res.ys2$M_star
cor.ys2_full.scaled <- 2*(cor.ys2_full-0.5)
# considering slope and time difference
cor.ys3_full <- w$w1*cor.S_star + w$w2*res.ys2$A_star2 + w$w3*res.ys2$M_star2
cor.ys3_full.scaled <- 2*(cor.ys3_full-0.5)

cor.yr1_full <- w$w1*cor.R_star + w$w2*res.ys1$A + w$w3*res.ys1$M
cor.yr1_full.scaled <- 2*(cor.yr1_full-0.5)
cor.yr2_full <- w$w1*cor.R_star + w$w2*res.ys2$A_star + w$w3*res.ys2$M_star
cor.yr2_full.scaled <- 2*(cor.yr2_full-0.5)
# considering slope and time difference
cor.yr3_full <- w$w1*cor.R_star + w$w2*res.ys2$A_star2 + w$w3*res.ys2$M_star2
cor.yr3_full.scaled <- 2*(cor.yr3_full-0.5)


# plot the important results

f_corrplot("cor.ys1", data=cor.ys1, order="original")
f_corrplot("cor.ys1", data=cor.ys1, order="hclust")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="original")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="hclust")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="original")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="hclust")

f_corrplot("cor.ys2", data=cor.ys2, order="original")
f_corrplot("cor.ys2", data=cor.ys2, order="hclust")

f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="original")
f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="hclust")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="original")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="hclust")


# scaling of the calculated correlation score
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="original")
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="hclust")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="original")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="hclust")

f_corrplot("cor.ys1_full.scaled", data=cor.ys1_full.scaled, order="original")
f_corrplot("cor.ys1_full.scaled", data=cor.ys1_full.scaled, order="hclust")
f_corrplot("cor.ys2_full.scaled", data=cor.ys2_full.scaled, order="original")
f_corrplot("cor.ys2_full.scaled", data=cor.ys2_full.scaled, order="hclust")
f_corrplot("cor.ys3_full.scaled", data=cor.ys3_full.scaled, order="original")
f_corrplot("cor.ys3_full.scaled", data=cor.ys3_full.scaled, order="hclust")

f_corrplot("cor.yr1_full.scaled", data=cor.yr1_full.scaled, order="original")
f_corrplot("cor.yr1_full.scaled", data=cor.yr1_full.scaled, order="hclust")
f_corrplot("cor.yr2_full.scaled", data=cor.yr2_full.scaled, order="original")
f_corrplot("cor.yr2_full.scaled", data=cor.yr2_full.scaled, order="hclust")
f_corrplot("cor.yr3_full.scaled", data=cor.yr3_full.scaled, order="original")
f_corrplot("cor.yr3_full.scaled", data=cor.yr3_full.scaled, order="hclust")


# extreme example where the difference between pearson and spearman matters
plot_cor_pair("Nos2", "Cxcl15")
plot_cor_pair("albumin", "Cyp2b10")

#---------------------------------------------
# Hierarchical clustering
#---------------------------------------------
# factor 1, factor 2, pearson, spearman, ys1_mean, ys1, ...
# -> plot the clusters and cluster members

# Perform hierarchical clustering on matrix
# The hclust function in R uses the complete linkage method for hierarchical clustering by default.
# This particular clustering method defines the cluster distance between two clusters to be the 
# maximum distance between their individual components.

# apply hirarchical clustering 
method = "yr3"
if (identical(method, "ys1")){
  cor.cluster <- cor.ys1_full.scaled  
}else if (identical(method, "ys2")){
  cor.cluster <- cor.ys2_full.scaled  
}else if (identical(method, "ys3")){
  cor.cluster <- cor.ys3_full.scaled  
}else if (identical(method, "yr1")){
  cor.cluster <- cor.yr1_full.scaled  
}else if (identical(method, "yr2")){
  cor.cluster <- cor.yr2_full.scaled  
}else if (identical(method, "yr3")){
  cor.cluster <- cor.yr3_full.scaled  
}

method = "yr1"
if (identical(method, "ys1")){
  cor.cluster <- cor.ys1_full.fil.scaled  
}else if (identical(method, "ys2")){
  cor.cluster <- cor.ys2_full.fil.scaled  
}else if (identical(method, "ys3")){
  cor.cluster <- cor.ys3_full.fil.scaled  
}else if (identical(method, "yr1")){
  cor.cluster <- cor.yr1_full.fil.scaled  
}else if (identical(method, "yr2")){
  cor.cluster <- cor.yr2_full.fil.scaled  
}else if (identical(method, "yr3")){
  cor.cluster <- cor.yr3_full.fil.scaled  
}


test <- dist(cor.cluster)

hc <- hclust(dist(cor.cluster)) 
corrplot(cor.cluster[hc$order, hc$order], order="original", method="square", type="full",tl.cex=0.3, tl.col="black")

# plot(hc)               # plot the dendrogram 
plot(hc, hang=-1)               # plot the dendrogram 

# cut tree into clusters
Ngroups = 8
rect.hclust(hc, k=Ngroups)
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)
groups.hc.order <- groups[hc$order]

dev.off()

# plot the groups
options$height = 3000
options$width = 3000

# Plot the clusters
for (k in 1:Ngroups){
  fname <- sprintf("%s_cluster_%s.png", method, k)
  png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  par(mfrow=c(N,N))
  for (name in names(g)){
    f_single_plot(name_A=name) 
  }
  par(mfrow=c(1,1))  
  dev.off()
}

install.packages('matrixStats')
library('matrixStats')

# mean plots for clusters 
# TODO: variance to the cluster
f_normalize_centering <- function(a){
  a.norm <- (a - mean(a))/(max(a, na.rm=TRUE) - min(a, na.rm=TRUE))
  return(a.norm)
}
options$height = 1600
options$width = 1600
fname <- sprintf("%s_cluster_overview.png", method)
png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
par(mfrow=c(ceiling(sqrt(Ngroups)),ceiling(sqrt(Ngroups))))
steps = 1:8
for (k in 1:Ngroups){
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  dgroup <- dmean[names(g)]

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


# better dendrogram
dend1 <- as.dendrogram(hc)
# Get the package:
install.packages("dendextend")
library(dendextend)

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
  