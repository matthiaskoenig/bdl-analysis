##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-21
#
#    This script contains the complete statistical analysis for the publication
#    Pathobiochemical signatures of cholestatic liver disease in bile duct ligated
#    mice (BMC Systems Biology).
#
#    The analysis consists of regression analysis of the various measured factors based
#    on different measures.
#    - Pearson Correlation 
#    - Spearman Correlation
#    - YS1 and YR1 (time course correlation)
#
##############################################################################################


# if filtering is selected the subset of factors being different 
# within the timecourse are selected
filter_by_anova = TRUE
if (filter_by_anova){
  # acceptance level
  p.accept = 0.05
  
  # how many rejected by adjusted p-value
  table(df.anova$p.holm>=p.accept) # 64 rejected / 90 accepted
  table(df.anova$p.value>=p.accept) # 19 rejected / 135 accepted
  
  # the accepted subset
  f.accept = df.anova$p.holm<p.accept
  
  # overwrite the full data set
  data.full <- data
  data <- data[, f.accept]
  factors.full <- factors
  factors <- factors[f.accept]
  dmean.full <- dmean
  dmean <- dmean[, f.accept]  
}



# plot of the data subset which is used for the correlation analysis
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")) 
heatmap.2(t(as.matrix(data)), col=col2(100), scale="row", Rowv=NULL, Colv=NULL,
          key=TRUE, trace="none", cexRow=0.5, keysize=0.8)

#---------------------------------------------
# Correlation analysis
#---------------------------------------------
# Calculation of standard correlation matrices based on Spearman
# and Pearson correlaton coefficients. 
# Hierarchical clustering is performed using complete linkage.
# The cor.pearson and cor.spearman are used as comparison models.

# Do the analysis either on the full dataset or the filtered subset




# correlation matrix
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")
cor.pearson.fil <- cor(data.fil, method="pearson", use="pairwise.complete.obs")
cor.spearman.fil <- cor(data.fil, method="spearman", use="pairwise.complete.obs")

install.packages("corrplot")
library(corrplot)
options <- list(width=2000, height=2000, res=200)

# helper function for correlation plot
f_corrplot <- function(name, data, order){
  fname <- sprintf("%s_%s.png", name, order)
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  png(filename=sprintf("../results/%s", fname), width=options$width, height=options$height, res=options$res)
  corrplot(data, order=order, hclust.method="complete", method="color", type="full", 
           tl.cex=0.3, tl.col="black", col=col2(10))
  dev.off()
}

# Spearman & Pearson correlation matrices
f_corrplot("cor.spearman", data=cor.spearman, order="original")
f_corrplot("cor.spearman", data=cor.spearman, order="hclust")
f_corrplot("cor.pearson", data=cor.pearson, order="original")
f_corrplot("cor.pearson", data=cor.pearson, order="hclust")
f_corrplot("cor.spearman.fil", data=cor.spearman.fil, order="original")
f_corrplot("cor.spearman.fil", data=cor.spearman.fil, order="hclust")
f_corrplot("cor.pearson.fil", data=cor.pearson.fil, order="original")
f_corrplot("cor.pearson.fil", data=cor.pearson.fil, order="hclust")



# --- YR1, YS1, YR2, YS2 correlation ----------
# Calculation of time-course based correlation measurements, namely ys1, ys2, yr1, yr2
# and the respective adaptions to the underlying datasets, i.e. using multiple
# repeat data with every time point measurment coming from an individual
# sample.
# In ys1, ys2, yr1, yr2 all calculations are performed on the mean time course data.
# The classical correlation components are replaced with the correlations calculated
# on the individual data points.

source("ys1_yr1.R")        # definition of ys1, ys2, yr1, yr2
source("ys1_yr1_tools.R")  # definition of ys1 and yr1

# calculation of ys1 and yr1 for mean data
w <- list(w1=0.5, w2=0.25, w3=0.25)
res.ys1 <- ys1.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1 <- res.ys1$value
cor.ys1.scaled <- 2*(cor.ys1-0.5)
res.ys1.fil <- ys1.df(dmean.fil, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1.fil <- res.ys1.fil$value
cor.ys1.fil.scaled <- 2*(cor.ys1.fil-0.5)

res.ys2 <- ys2.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2 <- res.ys2$value
cor.ys2.scaled <- 2*(cor.ys2-0.5)
res.ys2.fil <- ys2.df(dmean.fil, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2.fil <- res.ys2.fil$value
cor.ys2.fil.scaled <- 2*(cor.ys2.fil-0.5)


f_corrplot("cor.ys1", data=cor.ys1, order="original")
f_corrplot("cor.ys1", data=cor.ys1, order="hclust")
f_corrplot("cor.ys1.fil", data=cor.ys1.fil, order="original")
f_corrplot("cor.ys1.fil", data=cor.ys1.fil, order="hclust")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="original")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="hclust")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="original")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="hclust")

f_corrplot("cor.ys2", data=cor.ys2, order="original")
f_corrplot("cor.ys2", data=cor.ys2, order="hclust")
f_corrplot("cor.ys2.fil", data=cor.ys1.fil, order="original")
f_corrplot("cor.ys2.fil", data=cor.ys1.fil, order="hclust")

f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="original")
f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="hclust")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="original")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="hclust")

f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="original")
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="hclust")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="original")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="hclust")
f_corrplot("cor.ys1.fil.scaled", data=cor.ys1.fil.scaled, order="original")
f_corrplot("cor.ys1.fil.scaled", data=cor.ys1.fil.scaled, order="hclust")
f_corrplot("cor.ys2.fil.scaled", data=cor.ys2.fil.scaled, order="original")
f_corrplot("cor.ys2.fil.scaled", data=cor.ys2.fil.scaled, order="hclust")

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

# on filtered data
cor.ys1_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys1.fil$A + w$w3*res.ys1.fil$M
cor.ys1_full.fil.scaled <- 2*(cor.ys1_full.fil-0.5)
cor.ys2_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys2.fil$A_star + w$w3*res.ys2.fil$M_star
cor.ys2_full.fil.scaled <- 2*(cor.ys2_full.fil-0.5)
# considering slope and time difference
cor.ys3.fil <- w$w1*res.ys2.fil$S_star + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.ys3.fil.scaled <- 2*(cor.ys3.fil-0.5)
cor.ys3_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.ys3_full.fil.scaled <- 2*(cor.ys3_full.fil-0.5)

cor.yr1_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys1.fil$A + w$w3*res.ys1.fil$M
cor.yr1_full.fil.scaled <- 2*(cor.yr1_full.fil-0.5)

cor.yr2_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys2.fil$A_star + w$w3*res.ys2.fil$M_star
cor.yr2_full.fil.scaled <- 2*(cor.yr2_full.fil-0.5)
# considering slope and time difference
cor.yr3_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.yr3_full.fil.scaled <- 2*(cor.yr3_full.fil-0.5)

f_corrplot("cor.ys1_full.fil.scaled", data=cor.ys1_full.fil.scaled, order="hclust")
f_corrplot("cor.ys2_full.fil.scaled", data=cor.ys2_full.fil.scaled, order="hclust")
f_corrplot("cor.ys3_full.fil.scaled", data=cor.ys3_full.fil.scaled, order="hclust")

f_corrplot("cor.yr1_full.fil.scaled", data=cor.yr1_full.fil.scaled, order="hclust")
f_corrplot("cor.yr2_full.fil.scaled", data=cor.yr2_full.fil.scaled, order="hclust")
f_corrplot("cor.yr3_full.fil.scaled", data=cor.yr3_full.fil.scaled, order="hclust")


# name_A <- "Actb"
# name_B <- "Actb.x"
# extreme example where the difference between pearson and spearman matters
name_A <- "Nos2"
name_B <- "Cxcl15"
f_cor_pair_plot("Nos2", "Cxcl15")

f_cor_pair_plot("albumin", "Cyp2b10")

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
  