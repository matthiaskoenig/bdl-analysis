##############################################################################################
#
#    BDL - Decision Tree
#
#    Matthias Koenig
#    2015-07-20
#
#    Necessary to find relevant genes first, i.e. filtering.
#    
#    A one-way analyis of variance (ANOVA) was applied to filter genes showing
#    significant (padj<0.05) up- or down-regulation during the time course, using the Bonferonni
#    step-down procedure to correct for any artificial p-value inflation.
#
#   Andreas gap method:
#   Only those factors are considered where the values one time frame are
#   in a disjoint interval compared to the values of another time
#   frame. Between the intervals there is a gap, where the size of gap,
#   compared with the range of all values of both intervals --- the
#   relative gap --- measures the fitness of the factor to tell one time
#   point from another. If several factors are candidates the one with the
#   larger relative gap is perferred. The median of the gap is the
#   suggested splitting point between the two time points. This is called
#   gap method.
#
# Ordinal classification problem:
# R has (at least) several packages directed to ordinal regression. One of these is actually called Ordinal, but i haven't used it. I have used the Design Package in R for ordinal regression and i can certainly recommend it. Design contains a complete set of functions for solution, diagnostics, testing, and results presentation of ordinal regression problems via the Ordinal Logistic Model. Both Packages are available from CRAN) A step-by-step solution of an ordinal regression problem using the Design Package is presented on the UCLA Stats Site.
# http://www.ats.ucla.edu/stat/r/dae/ologit.htm
#
#
##############################################################################################

rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis")

# Load data for prediction
file.data <- file.path("..", "data", "bdl-data.Rdata")
load(file=file.data)

# The class to predict is the actual time point
rownames(data) <- paste(samples$time_fac, rownames(data), sep=" ")


# The possible predictors are all factors in the data set
# biocLite("ALL")
head(data)
# plot the predictors against classes
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  

library("RColorBrewer")
# display.brewer.all()


hc <- hclust(dist(cor.cluster)) 

# time colors
colorset <- brewer.pal(length(levels(samples$time_fac)), "Set2")
color.map <- function(time_f) {return(colorset[time_f])}
timeColors <- unlist(lapply(samples$time_fac, color.map))
# cluster colors
Ngroups = 8
# rect.hclust(hc, k=Ngroups)
groups <- cutree(hc, k=Ngroups)
colorset <- brewer.pal(Ngroups, "Set1")
color.map <- function(cluster_id) {return(colorset[cluster_id])}
clusterColors <- unlist(lapply(groups, color.map))


# TODO: use the cluster dendrogramm for reordering the factors
options <- list(width=1600, height=1600, res=200)
png(filename=sprintf("../results/%s", "factor_time_heatmap.png"), width=options$width, height=options$height, res=options$res)
heatmap.2(t(as.matrix(data)), col=col2(100), scale="row", Rowv=NULL, Colv=NULL,
          key=TRUE, trace="none", cexRow=0.5, keysize=0.8, ColSideColors=timeColors, RowSideColors=clusterColors)
dev.off()

png(filename=sprintf("../results/%s", "factor_time_heatmap_clusters.png"), width=options$width, height=options$height, res=options$res)
heatmap.2(t(as.matrix(data)), col=col2(100), scale="row", Rowv=as.dendrogram(hc), Colv=NULL, dendrogram="row",
          key=TRUE, trace="none", cexRow=0.5, keysize=0.8, ColSideColors=timeColors, RowSideColors=clusterColors)
dev.off()










heatmap.2(cor.cluster, col=col2(10), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5)


# The used predictors in the phase trees
# What was done in the original analysis
fac.tree <- c("CTGF", "alpha.SMA", "Tnfrsf1a", "Gstm1", "Il28b", "Fn1", "Il2")

options <- list(width=1600, height=1600, res=200)
png(filename="../results/tree_factors.png", width=options$width, height=options$height, res=options$res)
n.panel = ceiling(sqrt(length(fac.tree)))
par(mfrow=c(n.panel, n.panel))
for (name in fac.tree){
  f_single_plot(name)  
}
par(mfrow=c(1,1))
dev.off()

f_single_plot("Nr0b2")
sort(names(data))

library(rpart)
treedata <- data
treedata$class <- samples$time_fac
head(treedata)

# full formula
f = paste("class ~ ", paste(names(data), sep="", collapse=" + "), sep="")
# reudced formula to the factors in the decision tree so far
f = paste("class ~ ", paste(fac.tree, sep="", collapse=" + "), sep="")

# tree.fit <- rpart(formula=f, data=treedata, method="class", control=rpart.control(minsplit=1)) 
tree.fit <- rpart(formula=f, data=treedata, method="class", control=rpart.control(minsplit=5)) 
printcp(tree.fit)
print(tree.fit)
plot(tree.fit)
text(tree.fit)

#################################################



