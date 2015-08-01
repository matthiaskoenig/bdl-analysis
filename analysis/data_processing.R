require("calibrate")

#---------------------------------------------
# Read & preprocess data
#---------------------------------------------

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
rm(tmp, d, adme, antibodies, cytokines, fibrosis, histology, histology.processed)

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

# save processed data
file.data <- file.path("..", "data", "bdl-data.Rdata")
save(data, samples, dmean, dmean.time, file=file.data)
load(file=file.data)

#---------------------------------------------
# All data heatmap
#---------------------------------------------





#---------------------------------------------
# Single factor visualization
#---------------------------------------------
source("plots.R")

# Create plot of the first factor
plot_single_factor(1)

# Creates plots of all factors in data
plot_all_factors()


#--------------------------
# Actb controls
#--------------------------
# Actb was measured on all Fluidigm chips and serves as quality control of the 
# measurement/correlation analysis.
# Use pairwise correlation plots for control.
f_cor_pair_plot("Actb", "Por")

# Actb control figure
options <- list(width=1600, height=600, res=200)
png(filename="../results/Actb_control.png", width=options$width, height=options$height, res=options$res)
par(mfrow=c(1,3))
f_cor_pair_plot("Actb", "Actb.x", single_plots=FALSE)
f_cor_pair_plot("Actb", "Actb.y", single_plots=FALSE)
f_cor_pair_plot("Actb.x", "Actb.y", single_plots=FALSE)
par(mfrow=c(1,1))
dev.off()

# calculate the correlations
f_cor_pair_plot("Actb", "Actb.x")
f_cor_pair_plot("Actb", "Actb.y")
f_cor_pair_plot("Actb.x", "Actb.y")

cat('Actb Spearman : individual points\n')
actb.spearman <- cor(data.frame(Actb=data$Actb, 
                                Actb.x=data$Actb.x, 
                                Actb.y=data$Actb.y), method="spearman")
print(actb.spearman)
cat('Actb Spearman : mean\n')
actb.spearman.mean <- cor(data.frame(Actb=dmean$Actb, 
                                     Actb.x=dmean$Actb.x, 
                                     Actb.y=dmean$Actb.y), method="spearman")
print(actb.spearman.mean)
cat('Actb Pearson : individual points\n')
actb.pearson <- cor(data.frame(Actb=data$Actb, 
                               Actb.x=data$Actb.x, 
                               Actb.y=data$Actb.y), method="pearson")
print(actb.pearson)
cat('Actb Pearson : mean\n')
actb.pearson.mean <- cor(data.frame(Actb=dmean$Actb, 
                                    Actb.x=dmean$Actb.x, 
                                    Actb.y=dmean$Actb.y), method="pearson")
print(actb.pearson.mean)

f_cor_pair_plot("Actb", "Por")