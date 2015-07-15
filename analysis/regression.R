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

# Set sample ids as row numbers
rownames(data) <- data$sid

# create the ordered time factor
time <- data$time
time <- ordered(time, levels = c("0h", "6h", "12h", "18h", "30h", "2d", "5d", "14d"))
time

# remove non-factor columns which are not part of the correlation analysis
data <- subset(data, select = -c(sid,time) )
head(data)

#---------------------------------------------
# Mean factor timecourse
#---------------------------------------------
# mean based on times (aggregate)
library(reshape)
data2 <- data
data2$time <- samples$time
# rm NA in mean calculation
dmean <- aggregate(data2, list(data2$time), FUN=mean, na.rm=TRUE)
dmean.time <- dmean$time
dmean <- subset(dmean, select = -c(time, Group.1) )
rm(data2)

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
options = list(width=1600, height=800, res=150)

# factor names
factors <- colnames(data)

# Create plot of a single factor
# plot individual data points with number
plot_single_factor <- function(k){
  name <- factors[k]
  info <- get_probe_info(name)
  # print(name)
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
  points(dmean.time, dmean[,k], col="red", pch=15)
  lines(dmean.time, dmean[,k], col="red")
  
  # [B] plot as factor
  plot(time, data[,k], xlab="time", ylab=name, main=name, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(data[,k], na.rm=TRUE)*1.1))
  points(time, data[,k], col="black")
  points(time, data[,k], col=rgb(0,0,1,0.6), pch=16)
  
  points(1:nrow(dmean), dmean[,k], col="red", pch=15)
  lines(1:nrow(dmean), dmean[,k], col="red")

  par(mfrow=c(1,1))
  dev.off()
}
plot_single_factor(1)

# Create all the plots
k= 1
Nf = length(factors)
for (k in 1:Nf){
  cat(sprintf("%s / %s\n", k, Nf))
  plot_single_factor(k)
}

#---------------------------------------------
# Correlation analysis
#---------------------------------------------

# --- Pearson & Spearman correlation matrix -------------------------------------------------------------------
# calculation of correlation matrix
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")

# plotting of simple correlation matrices (original order) & based on hierarchical clustering.
# Using complete linkage method.

# corrplot options for clustering
# order = "hclust"
# hclust.method = "complete"
# corrplot(cor.spearman, method="square", type="full", 
#         tl.cex=0.3, tl.col="black",
#         order="hclust", hclust.method="complete", addrect=5)

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


# --- YR1 and YS1 correlation ------------------------------------------------------------------------------------

# Calculation is only possible on the mean time data
# => calculation has to be performed on the mean dataset
# cor.ys1 <- ys1.df(data, time, use="pairwise.complete.obs")



source("ys1_yr1.R")  # definition of ys1 and yr1
# calculation of ys1 and yr1 for mean data
res.ys1 <- ys1.df(data2.melt, data2.time, w1=0.25, w2=0.5, w3=0.25, use="pairwise.complete.obs")
cor.ys1 <- res.ys1$value

cor.yr1 <- yr1.df(data2.melt, data2.time, w1=0.25, w2=0.5, w3=0.25, use="pairwise.complete.obs")

f_corrplot("cor.ys1", data=cor.ys1, order="original")
f_corrplot("cor.ys1", data=cor.ys1, order="hclust")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="original")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="hclust")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="original")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="hclust")

which(res.ys1$A==1)
test <-res.ys1$A
levels(as.factor(test))
test[test<=0.6] = 0
corrplot(test, order="original", hclust.method="complete", method="square", type="full", 
         tl.cex=0.3, tl.col="black")



f_corrplot("cor.yr1", data=cor.yr1, order="original")
f_corrplot("cor.yr1", data=cor.yr1, order="hclust")

cor.ys1.scaled <- 2*(cor.ys1-0.5)
cor.yr1.scaled <- 2*(cor.yr1-0.5)
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="original")
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="hclust")
f_corrplot("cor.yr1.scaled", data=cor.yr1.scaled, order="original")
f_corrplot("cor.yr1.scaled", data=cor.yr1.scaled, order="hclust")


#--------------------------
# Actb controls
#--------------------------

install.packages("calibrate")
library(calibrate)
# use the textxy() function to add labels to the preexisting plot's points
# add labels for the total enrollment
  
# Actb control plots via pairwise correlation plots.
# Actb exists on all Fluidigm chips.


f_single_plot <- function(name_A){
  dA <- data[[name_A]]
  dA.mean <- dmean[[name_A]]
  plot(time, dA, xlab="time", ylab=name_A, main=name_A, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1))
  points(time, dA, col="black")
  points(time, dA, col=rgb(0,0,1,0.6), pch=16)
  
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
  dA.melt <- data2.melt[[name_A]]
  dB.melt <- data2.melt[[name_B]]
  
  value.max <- 1.05* max(max(dA), max(dB))
  plot(dA, dB, xlim=c(0, value.max), ylim=c(0, value.max), 
       main=sprintf("%s ~ %s", name_A, name_B),
       xlab=name_A, ylab=name_B)
  textxy(dA, dB, samples$time)
  points(dA.melt, dB.melt, pch=16, col=rgb(0,0,1,0.8) )
  lines(dA.melt, dB.melt, col=rgb(0,0,1,0.8) )
  textxy(dA.melt, dB.melt, data2.time, col="blue")
  abline(a=0, b=1, col="darkgray")
  
  # single plots
  if (single_plots){
    f_single_plot(name_A)
    f_single_plot(name_B)
  }
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
cor(data.frame(Actb=data2.melt$Actb, 
               Actb.x=data2.melt$Actb.x, 
               Actb.y=data2.melt$Actb.y), method="spearman")
cor(data.frame(Actb=data$Actb, 
               Actb.x=data$Actb.x, 
               Actb.y=data$Actb.y), method="pearson")
cor(data.frame(Actb=data2.melt$Actb, 
               Actb.x=data2.melt$Actb.x, 
               Actb.y=data2.melt$Actb.y), method="pearson")
#---------------------------------------------------------------------------


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



#---------------------------------------------
