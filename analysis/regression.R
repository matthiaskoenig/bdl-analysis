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
# install.packages("corrplot")
library(corrplot)

options$width=1600
options$height=1600
options$res=200

# Calculation of correlation scores
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")

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
