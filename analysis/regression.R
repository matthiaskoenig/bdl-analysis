##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-09
#
##############################################################################################

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
data <- subset(data, select = -c(sid,time) )
head(data)

#########################################################################
# Correlation analysis
#########################################################################
# First Correlogram Example
install.packages('corrgram')
library(corrgram)

corrgram(data, order=FALSE, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="BDL Correlation") 
