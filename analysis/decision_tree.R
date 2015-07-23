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
# sudo apt-get install r-cran-rgtk2
# install.packages("rpart.plot")
# install.packages("rattle")
# install.packages("party")

# require(party)
# require(rattle)

rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis")
source("plots.R")  # Load plot helpers

load(file=file.path("..", "data", "bdl-data.Rdata"))   # Load BDL data
load(file=file.path("..", "data", "bdl-anova.Rdata"))  # Load ANOVA results for filtering 

# The class to predict is the actual time point
rownames(data) <- paste(samples$time_fac, rownames(data), sep=" ")

# data subset with factors filtered via ANOVA
p.accepted <- df.anova$p.holm<0.05
data.fil <- data[, p.accepted]

# TODO: perform additional feature selection -> use clusters as features for the prediction


#--------------------------------------------------------------------------------------------

# Overview of predictors in the original analysis. 
# They were all good separators for certain time points and were selected by their
# maximal GAP.
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

#--------------------------------------------------------------------------------------------
# Standard decision trees 
require(rpart)
require(rpart.plot)

# add the class label to the predictors (time points/period are used as classes)
treedata <- data.fil
treedata$class <- samples$time_fac
head(treedata)

# create formula for tree classification
# full formula using all predictors of filtered data
formula.fil = paste("class ~ ", paste(colnames(data.fil), sep="", collapse=" + "), sep="")

# fit the tree with classification
# tree.fit <- rpart(formula=f, data=treedata, method="class", control=rpart.control(minsplit=1)) 
tree.fit <- rpart(formula=formula.fil, data=treedata, method="class", control=rpart.control(minsplit=5)) 
# print information of the tree fit
printcp(tree.fit)
print(tree.fit)
plot(tree.fit)
text(tree.fit)
prp(tree.fit, type=0, extra=101, yesno=TRUE)
# Problems with the divide and conquer algorithms. It maximizes the information gain on every split.


# try a regression on the times
treedata.reg <- data.fil
treedata.reg$regvalue <- log(samples$time+1)
head(treedata.reg)
formula.fil.reg = paste("regvalue ~ ", paste(colnames(data.fil), sep="", collapse=" + "), sep="")


sample(1:40, size=5)
tree.reg <- rpart(formula=formula.fil.reg, data=treedata.reg[sample(1:40, size=35), ], method="anova", control=rpart.control(minsplit=6, minbucket=2, cp=-1))
printcp(tree.reg)
print(tree.reg)
plot(tree.reg)
text(tree.reg)
prp(tree.reg, type=0, extra=101, yesno=TRUE)



# reudced formula to the factors in the decision tree so far
formula.old = paste("class ~ ", paste(fac.tree, sep="", collapse=" + "), sep="")
tree.old <- rpart(formula=formula.old, data=treedata, method="class", control=rpart.control(minsplit=3)) 
printcp(tree.old)
print(tree.old)
plot(tree.old)
text(tree.old)
prp(tree.old)

# rPartOrdinal -  An R package for deriving a classification tree for predicting an ordinal response
# Archer2010
# rpartOrdinal R package, which implements ordered twoing, the generalized Gini, and the ordinal impurity splitting
# methods. These splitting methods should be considered for use when deriving an ordinal
# response classification tree.

# install.packages("rpartOrdinal")
library(rpartOrdinal)
f = formula.fil
f
# 3.1. ordered twoing
otwoing.rpart <- rpart(f, data=treedata, method = twoing)
# otwoing.rpart <- rpart(f, data=treedata, method = twoing, control=rpart.control(minsplit=2))
printcp(otwoing.rpart)
plot(otwoing.rpart)
text(otwoing.rpart, pretty = TRUE)
prp(otwoing.rpart, type=0, extra=101, yesno=TRUE)
prp(tree.fit, type=0, extra=101, yesno=TRUE)

f_single_plot("Gstm1")
f_single_plot("Mki67")
f_single_plot("Cyp1a2")

# 3.2 ordinal impurity function
ordinal.rpart <- rpart(f, data = treedata, method = ordinal, control=rpart.control(minsplit=2))
plot(ordinal.rpart)
text(ordinal.rpart, pretty = TRUE)

# 3.3 Generalized Giny impurity
linear.loss.rpart <- rpart(f, method = "class", data=treedata, parms = list(loss=loss.matrix(method = "linear", treedata$class)),
                           control=rpart.control(minsplit=5, xval=40))
plot(linear.loss.rpart)
text(linear.loss.rpart, pretty = TRUE)

quad.loss.rpart <- rpart(f, method = "class", data=treedata[sample(1:40, size=35),], parms = list(loss=loss.matrix(method = "quad", treedata$class)),
                           control=rpart.control(minsplit=5, xval=40))
plot(quad.loss.rpart)
text(quad.loss.rpart, pretty = TRUE)

# 
# We note that another R package, party (Hothorn et al.
# 2009), can also be used to derive an ordinal conditional inference tree, where the variable
# selected for splitting a given node is determined using an inferential test (Hothorn et al. 2006). 
# These methods may prove useful when the dataset to be analyzed includes an ordinal
# response and the number of covariates exceeds the sample size.


# Better plots of the classification trees
# Plotting Classification Trees with the plot.rpart and rattle pckages
# install.packages("partykit")
# install.packages("party")
# install.packages("rpart.plot")
# sudo apt-get install r-cran-rgtk2
# install.packages("rattle")
install.packages("RColorBrewer")
library(rpart)				        # Popular decision tree algorithm
library(rattle)					# Fancy tree plot
library(rpart.plot)				# Enhanced tree plots
library(RColorBrewer)				# Color selection for fancy tree plot
library(party)					# Alternative decision tree algorithm
library(partykit)				# Convert rpart object to BinaryTree
library(rpart.plot)

tree.1 <- ordinal.rpart
# additional node labeling
#node.fun1 <- function(x, labs, digits, varlen)
#{
#  paste("t", x$frame$class)
#}
# prp(tree.1, type=0, extra=101, yesno=TRUE, node.fun=node.fun1)
prp(tree.1, type=0, extra=101, yesno=TRUE)
fancyRpartPlot(tree.1)
tree.1$frame


##########################

# The possible predictors are all factors in the data set
# biocLite("ALL")
head(data)
# plot the predictors against classes
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  

library("RColorBrewer")
# display.brewer.all()

# use the clustering information from the correlation analysis
hc <- hclust(dist(cor.cluster))

# time colors
colorset <- brewer.pal(length(levels(samples$time_fac)), "Set2")
color.map <- function(time_f) {return(colorset[time_f])}
timeColors <- unlist(lapply(samples$time_fac, color.map))
# cluster colors
Ngroups = 6
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

# rPartOrdinal -  An R package for deriving a classification tree for predicting an ordinal response
# Archer2010
install.packages("rpartOrdinal")

