##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-13
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

# TODO: better merge of data

#---------------------------------------------
# Read & preprocess data
#---------------------------------------------
rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis")

# Load the datasets
samples <- read.csv(file.path("..", "data", "01_samples.csv"), sep="\t")
histology <- read.csv(file.path("..", "data", "02_histology.csv"), sep="\t")

adme <- read.csv(file.path("..", "data", "03_fluidigm_ADME.csv"), sep="\t")
cytokines <- read.csv(file.path("..", "data", "04_fluidigm_cytokines.csv"), sep="\t")
fibrosis1 <- read.csv(file.path("..", "data", "05_fluidigm_fibrosis_01.csv"), sep="\t")
fibrosis2 <- read.csv(file.path("..", "data", "06_fluidigm_fibrosis_02.csv"), sep="\t")

antibodies <- read.csv(file.path("..", "data", "07_antibodies.csv"), sep="\t")

# For the fibrosis panel two chips were measured. The mean value of the two chips is
# used for analysis.
fibrosis <- (fibrosis1 + fibrosis2) / 2
fibrosis$time <- fibrosis1$time
rm(fibrosis1, fibrosis2)

# prepare histology & antibody data
d <- list()
d$gldh <- data.frame(histology[, c("sid_GLDH", "GLDH")])
names(d$gldh) <- c("sid", "GLDH")
d$alt <- data.frame(histology[, c("sid_ALT", "ALT")])
names(d$alt) <- c("sid", "ALT")
d$bilirubin <- data.frame(histology[, c("sid_Bilirubin", "Bilirubin")])
names(d$bilirubin) <- c("sid", "bilirubin")
d$albumin <- data.frame(histology[, c("sid_albumin", "albumin")])
names(d$albumin) <- c("sid", "albumin")
d$hc <- data.frame(histology[, c("sid_BrdU_HC", "BrdU_HC")])
names(d$hc) <- c("sid", "BrdU_HC")
d$nhc <- data.frame(histology[, c("sid_BrdU_NHC", "BrdU_NHC")])
names(d$nhc) <- c("sid", "BrdU_NHC")
d$kupffer <- data.frame(histology[, c("sid_BrdU_Kupffer", "BrdU_Kupffer")])
names(d$kupffer) <- c("sid", "BrdU_Kupffer")
d$hsc <- data.frame(histology[, c("sid_BrdU_HSC", "BrdU_HSC")])
names(d$hsc) <- c("sid", "BrdU_HSC")
d$siriusRed <- data.frame(histology[, c("sid_BrdU_SiriusRed", "BrdU_SiriusRed")])
names(d$siriusRed) <- c("sid", "BrdU_SirirusRed")
d$bileInfarcts <- data.frame(histology[, c("sid_bileInfarcts", "bileInfarcts")])
names(d$bileInfarcts) <- c("sid", "bileInfarcts")
summary(d)

# Merge the histological & antibody datasets on sample ids
tmp <- merge(d$gldh, d$alt, by="sid")
tmp <- merge(tmp, d$bilirubin, by="sid")
tmp <- merge(tmp, d$albumin, by="sid")
tmp <- merge(tmp, d$hc, by="sid")
tmp <- merge(tmp, d$nhc, by="sid")
tmp <- merge(tmp, d$kupffer, by="sid")
tmp <- merge(tmp, d$hsc, by="sid")
tmp <- merge(tmp, d$siriusRed, by="sid")
tmp <- merge(tmp, d$bileInfarcts, by="sid")
histology.processed <- tmp
head(histology.processed)

# remove pid from antibodies
antibodies <- subset(antibodies, select = -c(pid) )
head(antibodies)

# Merge histological & antibody data with the fluidigm data
tmp <- merge(adme, cytokines, by = c('sid', 'time'))
tmp <- merge(tmp, fibrosis, by = c('sid', 'time'))
tmp <- merge(tmp, histology.processed, by=c('sid'))
tmp <- merge(tmp, antibodies, by=c('sid'))
data <- tmp

# Set sample ids as row numbers for data and samples
rownames(data) <- data$sid
rownames(samples) <- samples$sid

# create the ordered time factor in the samples
samples$time_fac <- ordered(data$time, levels = c("0h", "6h", "12h", "18h", "30h", "2d", "5d", "14d"))
samples$time_point <- as.vector(t(matrix(rep(seq(1,8), 5), nrow=8, ncol=5)))
# create the repeat info
samples$repeats <- (seq(from=0, to=(nrow(samples)-1))%%5) +1

# remove non-factor columns which are not part of the correlation analysis
data <- subset(data, select = -c(sid,time) )
rm(tmp)

# factor names
factors <- colnames(data)

head(data)
head(samples)

#---------------------------------------------
# Data reshaping
#---------------------------------------------
# data reformating for correlation algorithms.

# BDL mean data via aggregation on times
bdl_mean_data <- function(data){  
  library(reshape)
  data2 <- data
  data2$time <- samples$time
  # rm NA in mean calculation
  dmean <- aggregate(data2, list(data2$time), FUN=mean, na.rm=TRUE)
  dmean.time <- dmean$time
  dmean <- subset(dmean, select = -c(time, Group.1))
  return( list(dmean=dmean, dmean.time=dmean.time) )
}

# Mean of factors for the individual time points
tmp <- bdl_mean_data(data)
dmean <- tmp$dmean
dmean.time <- tmp$dmean.time
rm(tmp)

head(dmean)
head(dmean.time)

# List of data matrices
bdl_data_matrices <- function(data, time_pts, Nrepeats=5){
  data_list <- list()
  for (name in names(data)){
    # important to fill in the right order !
    data_list[[name]] <- matrix(data[[name]], nrow=length(dmean.time), ncol=Nrepeats, byrow=TRUE)
    colnames(data_list[[name]]) <- paste('R', 1:5, sep="")
    rownames(data_list[[name]]) <- time_pts
  }
  names(data_list) <- names(data)
  return(data_list)
}
data_list <- bdl_data_matrices(data, time_pts=levels(time))
summary(data_list)
data_list[1]


#---------------------------------------------
# Fluidigm annotation information
#---------------------------------------------
# helper function with information for fluidigm probes
probes <- read.csv(file.path("..", "data", "probe_mapping.csv"), stringsAsFactors=FALSE, sep="\t")
head(probes)

get_probe_info <- function(gene_id){
  info <- list()
  idx <- which(probes$Gene==gene_id)
  if (!is.null(idx) & length(idx)>0){
    info$Protein.name <- probes$Protein.names[idx[1]]
    info$Chip <- probes$Chip[idx[1]]
    info$Entry <- probes$Entry[idx[1]]
    info$Gene.name <- probes$Gene.names[idx[1]]
  } else {
    info$Protein.name <- NULL
  }
  return(info)
}
get_probe_info('Ppara')

#---------------------------------------------
# Single factor analysis
#---------------------------------------------
# plot options
options = list(width=1600, height=800, res=150)

# Create plot of a single factor
# plot individual data points with number
plot_single_factor <- function(k){
  name <- factors[k]
  dA <- data[,k]  # factor data
  info <- get_probe_info(name)  # fluidigm probe annotation (if existing)
  
  fname <- paste("../results/factors/", sprintf("%03d", k), "_", name, ".png", sep="")
  png(filename=fname, width=options$width, height=options$height, res=options$res)
  par(mfrow=c(1,2))
  
  # [A] plot with time
  plot(samples$time_fac, dA, at=sort(as.numeric(levels(as.factor(samples$time)))), col="blue",
       xlab="time [h]", ylab=name, main=name,
       ylim=c(0, max(dA, na.rm=TRUE)*1.1 ))
  
  points(samples$time, dA, col="black")
  points(samples$time, dA, col=rgb(0,0,1,0.6), pch=16)
  if (!is.null(info$Protein.name)){
    text(x=140, y=max(dA, na.rm=TRUE)*1.08, 
         labels=info$Protein.name, cex=0.8)
  }
  points(dmean.time, dmean[,k], col="red", pch=15)
  lines(dmean.time, dmean[,k], col="red")
  
  # [B] plot as factor
  plot(samples$time_fac, dA, xlab="time", ylab=name, main=name, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(data[,k], na.rm=TRUE)*1.1))
  points(samples$time_fac, data[,k], col="black")
  points(samples$time_fac, data[,k], col=rgb(0,0,1,0.6), pch=16)
  
  points(1:nrow(dmean), dmean[,k], col="red", pch=15)
  lines(1:nrow(dmean), dmean[,k], col="red")

  par(mfrow=c(1,1))
  dev.off()
}
plot_single_factor(1)

# Create all single factor plots
plot_all_factors <- function(){
  Nf = length(factors)
  for (k in 1:Nf){
    cat(sprintf("%s / %s\n", k, Nf))
    plot_single_factor(k)
  }  
}
plot_all_factors()


#---------------------------------------------
# Correlation analysis
#---------------------------------------------

# --- Pearson & Spearman correlation matrix -------------------------------------------------------------------
# calculation of correlation matrix
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")

# plotting of simple correlation matrices (original order) & based on hierarchical clustering.
# Using complete linkage method.

library(corrplot)
options <- list(width=1600, height=1600, res=200)

# helper function for correlation plot
f_corrplot <- function(name, data, order){
  fname <- sprintf("%s_%s.png", name, order)
  png(filename=sprintf("../results/%s", fname), width=options$width, height=options$height, res=options$res)
  corrplot(data, order=order, hclust.method="complete", method="square", type="full", 
           tl.cex=0.3, tl.col="black")
  dev.off()
}
# Spearman & Pearson plots
f_corrplot("cor.spearman", data=cor.spearman, order="original")
f_corrplot("cor.spearman", data=cor.spearman, order="hclust")
f_corrplot("cor.pearson", data=cor.pearson, order="original")
f_corrplot("cor.pearson", data=cor.pearson, order="hclust")


# --- YR1, YS1, YR2, YS2 correlation ------------------------------------------------------------------------------------

# Calculation is performed on the mean timecourse dataset
source("ys1_yr1.R")  # definition of ys1 and yr1
source("ys1_yr1_tools.R")  # definition of ys1 and yr1

# calculation of ys1 and yr1 for mean data
w <- list(w1=0.5, w2=0.25, w3=0.25)
res.ys1 <- ys1.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1 <- res.ys1$value
res.ys2 <- ys2.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2 <- res.ys2$value

res.ys2$A_star2

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

cor.ys1.scaled <- 2*(cor.ys1-0.5)
cor.ys2.scaled <- 2*(cor.ys2-0.5)

f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="original")
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="hclust")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="original")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="hclust")

#--------------------------
# Actb controls
#--------------------------

# install.packages("calibrate")
library(calibrate)
# use the textxy() function to add labels to the preexisting plot's points
# add labels for the total enrollment
  
# Actb control plots via pairwise correlation plots.
# Actb exists on all Fluidigm chips.

f_single_plot <- function(name_A){
  dA <- data[[name_A]]
  dA.mean <- dmean[[name_A]]
  plot(samples$time_fac, dA, xlab="time", ylab=name_A, main=name_A, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1))
  points(samples$time_fac, dA, col="black")
  points(samples$time_fac, dA, col=rgb(0,0,1,0.6), pch=16)
  # plot the repeat number
  dA.text <- dA
  dA[is.na(dA)] <- -1
  textxy(samples$time_point, dA, samples$repeats, col="black", cex=1.1)
  
  points(1:nrow(dmean), dA.mean, col="red", pch=15)
  lines(1:nrow(dmean), dA.mean, col="red")
}

f_cor_pair_plot <- function(name_A, name_B, single_plots=TRUE){
  if (single_plots){
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
  }
  # correlation plot
  dA <- data[[name_A]]
  dB <- data[[name_B]]
  # dA <- dA[!is.na(dA)]
  # d <- dA[!is.na(dA)]
  dA.mean <- dmean[[name_A]]
  dB.mean <- dmean[[name_B]]
  
  value.max <- 1.05* max(max(dA), max(dB))
  plot(dA, dB, xlim=c(0, value.max), ylim=c(0, value.max), 
       main=sprintf("%s ~ %s", name_A, name_B),
       xlab=name_A, ylab=name_B)
  textxy(dA, dB, samples$time)
  points(dA.mean, dB.mean, pch=16, col=rgb(0,0,1,0.8) )
  lines(dA.mean, dB.mean, col=rgb(0,0,1,0.8) )
  textxy(dA.mean, dB.mean, dmean.time, col="blue")
  abline(a=0, b=1, col="darkgray")
  
  # single plots
  if (single_plots){
    f_single_plot(name_A)
    f_single_plot(name_B)
  }
  layout(matrix(c(1), 1, 1, byrow = TRUE)) 
}

options <- list(width=1600, height=600, res=200)
png(filename="../results/Actb_control.png", width=options$width, height=options$height, res=options$res)
par(mfrow=c(1,3))
f_cor_pair_plot("Actb", "Actb.x", single_plots=FALSE)
f_cor_pair_plot("Actb", "Actb.y", single_plots=FALSE)
f_cor_pair_plot("Actb.x", "Actb.y", single_plots=FALSE)
par(mfrow=c(1,1))
dev.off()

# calculate the correlations
cor(data.frame(Actb=data$Actb, 
               Actb.x=data$Actb.x, 
               Actb.y=data$Actb.y), method="spearman")
cor(data.frame(Actb=dmean$Actb, 
               Actb.x=dmean$Actb.x, 
               Actb.y=dmean$Actb.y), method="spearman")
cor(data.frame(Actb=data$Actb, 
               Actb.x=data$Actb.x, 
               Actb.y=data$Actb.y), method="pearson")
cor(data.frame(Actb=dmean$Actb, 
               Actb.x=dmean$Actb.x, 
               Actb.y=dmean$Actb.y), method="pearson")
#---------------------------------------------------------------------------

source("ys1_yr1.R")  # definition of ys1 and yr1
source("ys1_yr1_tools.R")  # definition of ys1 and yr1
# Full pearson correlation with mean slope & min/max component
cor.s <- ( cor(data, method="spearman", use="pairwise.complete.obs") + 1 )/2
cor.r <- ( cor(data, method="pearson", use="pairwise.complete.obs") + 1 )/2

cor.test <- 0.5*cor.s + 0.25*res.ys2$A_star + 0.25*res.ys2$M_star
cor.test.scaled <- 2*(cor.test-0.5)

cor.test.r <- 0.4*cor.r + 0.4*res.ys2$A_star + 0.2*res.ys2$M_star
cor.test.r.scaled <- 2*(cor.test.r-0.5)

cor.test.r2 <- 0.4*cor.r + 0.4*res.ys2$A_star2 + 0.2*res.ys2$M_star
cor.test.r2.scaled <- 2*(cor.test.r2-0.5)

summary(cor.test.r2.scaled <- 2*(cor.test.r2-0.5))
res.ys2$M_star2

son.fM_star(dmean[["Actb"]], dmean[["Actb.x"]], dmean.time)
son.fM_star2(dmean[["Actb"]], dmean[["Actb.x"]], dmean.time)
f_cor_pair_plot("Actb", "Actb.x")

summary(dmean["Actb.x"])


# cor.test.r <- 0.5*r.test + 0.25*res.ys2$A_star + 0.25*res.ys2$M_star
# cor.test.r.scaled <- 2*(cor.test.r-0.5)



f_corrplot("cor.test.scaled", data=cor.test.scaled, order="original")
f_corrplot("cor.test.scaled", data=cor.test.scaled, order="hclust")

f_corrplot("cor.test.r.scaled", data=cor.test.r.scaled, order="original")
f_corrplot("cor.test.r.scaled", data=cor.test.r.scaled, order="hclust")



# name_A <- "Actb"
# name_B <- "Actb.x"
# extreme example where 
name_A <- "Nos2"
name_B <- "Cxcl15"

f_cor_pair_plot(name_A, name_B)
# ys1(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=w$w1, w2=w$w2, w3=w$w3)
ys2(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=w$w1, w2=w$w2, w3=w$w3)
cor(x=data[[name_A]], data[[name_B]], method="spearman")
cor(x=dmean[[name_A]], dmean[[name_B]], method="spearman")



son.fS_star(a=data[[name_A]], b=data[[name_B]])


#---------------------------------------------
# Hierarchical clustering
#---------------------------------------------
# factor 1, factor 2, pearson, spearman, ys1_mean, ys1, ...
# -> plot the clusters and cluster members

# Perform hierarchical clustering on matrix
# The hclust function in R uses the complete linkage method for hierarchical clustering by default.
# This particular clustering method defines the cluster distance between two clusters to be the 
# maximum distance between their individual components.

# hc <- hclust(dist(cor.ys2))  # apply hirarchical clustering 
# hc <- hclust(dist(cor.test.scaled))  # apply hirarchical clustering 
hc <- hclust(dist(cor.test.r2.scaled))  # apply hirarchical clustering 

corrplot(cor.test.r2.scaled[hc$order, hc$order], order="original", method="square", type="full",tl.cex=0.3, tl.col="black")

# plot(hc)               # plot the dendrogram 
plot(hc, hang=-1)               # plot the dendrogram 

# cut tree into clusters
Ngroups = 8
rect.hclust(hc, k=Ngroups)
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)

groups.hc.order <- groups[hc$order]
groups.hc.order

# plot the groups
options$height = 3000
options$width = 3000

for (k in 1:Ngroups){
  fname <- sprintf("%s_cluster_%s.png", "ys2", k)
  png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  print(N)
  par(mfrow=c(N,N))
  for (name in names(g)){
    f_single_plot(name_A=name) 
  }
  par(mfrow=c(1,1))  
  dev.off()
}

f_cor_pair_plot("Nos2", "Hk2")

# make the mean plots (TODO: add the mean and variance of the cluster)

par(mfrow=c(2,4))
for (k in 1:Ngroups){
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))

  plot(1, type="n", xlab="", ylab="", xlim=c(1, 8), ylim=c(-1, 1), main=sprintf("Cluster %s", k))
  for (name in names(g)){
    dA = dmean[[name]]
    points(1:8, (dA - mean(dA))/(max(dA)-min(dA)), pch=16, col="black")
    lines(1:8, (dA - mean(dA))/(max(dA)-min(dA)), col="blue")
  }
}
par(mfrow=c(1,1))



# TODO: create the list of all pairwise correlations

  