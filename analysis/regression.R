##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-13
#
##############################################################################################

#---------------------------------------------
# Read & preprocess data
#---------------------------------------------
# There are two chips measured for the fibrosis score. 
# In the first version of this analysis only the first panel is used.
# TODO: necessary to include both panels for analysis
rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis")

# Load the datasets
samples <- read.csv(file.path("..", "data", "01_samples.csv"), sep="\t")
histology <- read.csv(file.path("..", "data", "02_histology.csv"), sep="\t")

adme <- read.csv(file.path("..", "data", "03_fluidigm_ADME.csv"), sep="\t")
cytokines <- read.csv(file.path("..", "data", "04_fluidigm_cytokines.csv"), sep="\t")
fibrosis1 <- read.csv(file.path("..", "data", "05_fluidigm_fibrosis_01.csv"), sep="\t")
fibrosis2 <- read.csv(file.path("..", "data", "06_fluidigm_fibrosis_02.csv"), sep="\t")

# prepare the histology data
names(histology)
gldh <- data.frame(histology[, c("sid_GLDH", "GLDH")])
names(gldh) <- c("sid", "GLDH")
alt <- data.frame(histology[, c("sid_ALT", "ALT")])
names(alt) <- c("sid", "ALT")
bilirubin <- data.frame(histology[, c("sid_Bilirubin", "Bilirubin")])
names(bilirubin) <- c("sid", "bilirubin")
albumin <- data.frame(histology[, c("sid_albumin", "albumin")])
names(albumin) <- c("sid", "albumin")
hc <- data.frame(histology[, c("sid_BrdU_HC", "BrdU_HC")])
names(hc) <- c("sid", "BrdU_HC")
nhc <- data.frame(histology[, c("sid_BrdU_NHC", "BrdU_NHC")])
names(nhc) <- c("sid", "BrdU_NHC")
kupffer <- data.frame(histology[, c("sid_BrdU_Kupffer", "BrdU_Kupffer")])
names(kupffer) <- c("sid", "BrdU_Kupffer")
hsc <- data.frame(histology[, c("sid_BrdU_HSC", "BrdU_HSC")])
names(hsc) <- c("sid", "BrdU_HSC")
siriusRed <- data.frame(histology[, c("sid_BrdU_SiriusRed", "BrdU_SiriusRed")])
names(siriusRed) <- c("sid", "BrdU_SirirusRed")

tmp <- merge(gldh, alt, by="sid")
tmp <- merge(tmp, bilirubin, by="sid")
tmp <- merge(tmp, albumin, by="sid")
tmp <- merge(tmp, hc, by="sid")
tmp <- merge(tmp, nhc, by="sid")
tmp <- merge(tmp, kupffer, by="sid")
tmp <- merge(tmp, hsc, by="sid")
histology.processed <- merge(tmp, siriusRed, by="sid")
head(histology.processed)

# Merge datasets on sid, time
tmp <- merge(adme, cytokines, by = c('sid', 'time'))
tmp <- merge(tmp, fibrosis1, by = c('sid', 'time'))
# currently only the first data set used
# TODO: method to handle both fibrosis data sets (mean?)
# tmp <- merge(tmp, fibrosis2, by = c('sid', 'time'))

data <- merge(tmp, histology.processed, by=c('sid'))
names(data)

# remove columns not in the correlation analysis
rownames(data) <- data$sid
head(data)

# the ordered time factors
time <- data$time
time <- ordered(time, levels = c("0h", "6h", "12h", "18h", "30h", "2d", "5d", "14d"))
time

data <- subset(data, select = -c(sid,time) )
head(data)

#---------------------------------------------
# Fluidigm annotation information
#---------------------------------------------
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
options = list()
options$width=1600
options$height=800
options$res=150

# factor names
factors <- colnames(data)

# Create plot of a single factor
# plot individual data points with number
# TODO: add estimated curves (splines)
plot_single_factor <- function(k){
  name <- factors[k]
  info <- get_probe_info(name)
  fname <- paste("../results/factors/", sprintf("%03d", k), "_", name, ".png", sep="")
  png(filename=fname, width=options$width, height=options$height, res=options$res)
  
  par(mfrow=c(1,2))
  
  # [A] plot with time
  plot(time, data[,k], at=sort(as.numeric(levels(as.factor(samples$time)))), col="blue",
       xlab="time [h]", ylab=name, main=name,
       ylim=c(0, max(data[,k], na.rm=TRUE)*1.1 ))
  
  points(samples$time, data[,k], col="black")
  points(samples$time, data[,k], col=rgb(0,0,1,0.6), pch=16)
  if (!is.null(info$Protein.name)){
    text(x=140, y=max(data[,k], na.rm=TRUE)*1.08, 
         labels=info$Protein.name, cex=0.8)
  }
  
  # [B] plot as factor
  plot(time, data[,k], xlab="time", ylab=name, main=name, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(data[,k], na.rm=TRUE)*1.1))
  points(time, data[,k], col="black")
  points(time, data[,k], col=rgb(0,0,1,0.6), pch=16)

  par(mfrow=c(1,1))
  dev.off()
}
# plot_single_factor(1)

# Create all the plots
for (k in 1:length(factors)){
  print(k)
  plot_single_factor(k)
}

#---------------------------------------------
# Correlation analysis
#---------------------------------------------

# ys1 and yr1 definition [Son2007]
# "A modfied correlation coefficient based similarity measure for clustering
#   time-course gene expression data".

# Slope between two adjacent measurments
slope <- function(x1,x2,t1,t2){
  s = ((x2 - x1) / (t2 - t1))  
}

# Sign of slope
ys1_L <- function(x1,x2,t1,t2)
  L = sign(slope(x1,x2,t1,t2));
end


# Relative count of the identical slopes.
ys1_fA <- function(a,b,time_pts){
  N = length(a)
  count = 0
  for (i in 1:(N-1)){
    
    if (ys1_L(a[i],a[i+1],time_pts[i],time_pts[i+1]) == ys1_L(b[i],b[i+1],time_pts[i],time_pts[i+1])){
       count = count + 1
    }
  }
  r = count / (N-1)
}

# Adapted Spearman correlation S*
ys1_fS_star <- function(a,b){
  r = (cor(x=x1,y=x2, method='spearman') + 1) / 2  
}

# Compare the timepoint of minimal and maximal value.
# ? how to handle NA -> modify to use the pairwise complete observations
ys1_fM = function(a,b){
  # find index of max and min
  idx_max_a = which.max(a)
  idx_max_b = which.max(b)
  idx_min_a = which.min(a)
  idx_min_b = which.min(b)
  
  # compare max and min indices
  if ( (idx_max_a == idx_max_b) && (idx_min_a == idx_min_b) ){
     r = 1.0
  } else if ( (idx_max_a == idx_max_b) || (idx_min_a == idx_min_b) ){
     r = 0.5
  } else if ( (idx_max_a != idx_max_b) && (idx_min_a != idx_min_b) ){
     r = 0.0
  }
  return(r)
}

ys1 = function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  # handle the NAs
  if (identical(use, "pairwise.complete.obs")){
    df <- data.frame(a, b, time_pts)
    df <- df[complete.cases(df),]
    a <- df$a
    b <- df$b
    time_pts <- df$time_pts
  }  
  r = (w1 * ys1_fS_star(a,b)) + (w2 * ys1_fA(a,b,time_pts)) + (w3 * ys1_fM(a,b))
  return(r)
}

# test data for implementation
# yR1(x1, x2)
time_pts <- seq(0,4)
x1 <- c(2,3,6,4,7)
x2 <- c(1,2,3,5,3)
y1 <- c(4,3,6,2,7)
y2 <- c(5,2,3,1,3)
par(mfrow=c(1,2))
plot(time_pts, x1, col="blue", pch=16, type="o", ylim=c(0, max(max(x1), max(x2))),
     main = "x1, x2", xlab="values", ylab="time")
points(time_pts, x2, col="red", pch=15, type="o")
plot(time_pts, y1, col="blue", pch=16, type="o", ylim=c(0, max(max(y1), max(y2))),
     main = "y1, y2", xlab="values", ylab="time")
points(time_pts, y2, col="red", pch=15, type="o")
par(mfrow=c(1,1))

# simple correlation values
cor(x=x1,y=x2, method='spearman') # 0.667
cor(x=x1,y=x2, method='pearson') # 0.439

# YS1
ys1(x1, x2, time_pts)
ys1(x1, x2, time_pts) == 0.583
ys1(y1, y2, time_pts)
ys1(y1, y2, time_pts) == 0.833

# TODO: YR1 values
# yr1(x1, x2, time_pts)
# yr1(x1, x2, time_pts) == 0.555
# yr1(y1, y2, time_pts)
# yr1(y1, y2, time_pts) == 0.805

# Calculate ys1 for a given data frame
ys1.df <- function(data, time_pts, use="pairwise.complete.obs"){
  N <- ncol(data)
  cor.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(cor.mat) <- names(data)
  rownames(cor.mat) <- names(data)
  for (k in 1:N){
    for (i in 1:N){
       cor.mat[k,i] = ys1(data[,k], data[,i], time_pts, use=use)
    }
  }
  return(cor.mat)
}

# TODO: this is not possible on the complete data set, 
# => calculation has to be performed on the mean dataset
cor.ys1 <- ys1.df(data, time, use="pairwise.complete.obs")

# mean based on times (aggregate)
library(reshape)
data2 <- data
data2$time <- samples$time

data2.melt <- aggregate(data2, list(data2$time), FUN=mean)
data2.time <- data2.melt$time
data2.melt <- subset(data2.melt, select = -c(time, Group.1) )
cor.ys1 <- ys1.df(data2.melt, data2.time, use="pairwise.complete.obs")

# install.packages("corrplot")
options$width=1600
options$height=1600
options$res=200
library(corrplot)

cor.ys1.scaled <- 2*(cor.ys1-0.5)

# TODO: bug -> check the act results !!!! Either naming or calculation of the correlation
print('----')
ys1(data2.melt$Actb, data2.melt$Act.B.x, data2.time, use="pairwise.complete.obs")
ys1(data2.melt$Actb, data2.melt$Act.B.y, data2.time, use="pairwise.complete.obs")
ys1(data2.melt$Act.B.x, data2.melt$Act.B.y, data2.time, use="pairwise.complete.obs")
# ??? does not depend at all on second argument ??? problem


png(filename="../results/cor.ys1-scaled_original.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.ys1.scaled, order="original", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()
png(filename="../results/cor.ys1-scaled_hclust.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.ys1.scaled, order="hclust", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()





corrplot(cor.spearman, order="hclust", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()





# Calculation of correlation scores
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")

data[,143]
names(data)[143]
ys1(a=data[, 1], b=data[,143], time_pts=samples$time, use="pairwise.complete.obs")

# Create the




# Spearman correlation plots
png(filename="../results/cor.spearman_original.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.spearman, order="original", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()

png(filename="../results/cor.spearman_hclust.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.spearman, order="hclust", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()

# Pearson
png(filename="../results/cor.pearson_original.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.pearson, order="original", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()

png(filename="../results/cor.pearson_hclust.png", width=options$width, height=options$height, 
    res=options$res)
corrplot(cor.pearson, order="hclust", method="square", type="full", 
         tl.cex=0.3, tl.col="black", # label settings
)
dev.off()

#---------------------------------------------





# First Correlogram Example
# install.packages('corrgram')
# install.packages('lattice')
# install.packages('ellipse')

library(corrgram)
library(lattice)
library(ellipse)

colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
colors <- colorfun(seq(from=-0.5, to=0.5, length.out=100))

# plotcorr(cor.pearson, col=rgb(colorfun((cor.pearson+1)/2), maxColorValue=255),
#          mar = c(0.1, 0.1, 0.1, 0.1))

levelplot(cor.spearman, main="correlation matrix", xlab="", ylab="")

rgb.palette <- colorRampPalette(c("red", "white", "green"), space = "rgb")
levelplot(cor.pearson, xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))

# TODO: [-1, 1] correlation, left upper corner starting
rgb.palette <- colorRampPalette(c("green", "white", "red"), space = "rgb")
levelplot(cor.spearman, xlab="", ylab="", scales=list(x=list(rot=90)), col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))






# plot all single timecourses



library(ggplot2)
library(plyr)
library(reshape)

cor.pearson.m <- melt(cor.pearson)
ggplot(z.m, aes(X1, X2, fill = value)) + geom_tile() + 
  scale_fill_gradient2(low = "blue",  high = "yellow")
