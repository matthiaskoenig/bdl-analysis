#-------------------------------------------------------------------------------
# Read & preprocess data
#
# The individual raw data files are loaded and combined in a single data frame.
# Factor and timecourse information is added.
# All processed datasets are stored as Rdata files in combination with the 
# packages.
#
# Processed data is made available in the package via
#     library(BDLanalysis)
#     data(BDLdata)
#     data(BDLsamples)
#-------------------------------------------------------------------------------
rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/BDLanalysis")

samples <- read.csv('inst/extdata/01_samples.csv', sep="\t")
histology <- read.csv('inst/extdata/02_histology.csv', sep="\t")
adme <- read.csv('inst/extdata/03_fluidigm_ADME.csv', sep="\t")
cytokines <- read.csv('inst/extdata/04_fluidigm_cytokines.csv', sep="\t")
fibrosis1 <- read.csv('inst/extdata/05_fluidigm_fibrosis_01.csv', sep="\t")
fibrosis2 <- read.csv('inst/extdata/06_fluidigm_fibrosis_02.csv', sep="\t")
antibodies <- read.csv('inst/extdata/07_antibodies.csv', sep="\t")

# Two repeats of the fibrosis fluidigm chip were measured. 
# The mean value of the two chips is used for analysis.
fibrosis <- (fibrosis1 + fibrosis2) / 2
fibrosis$time <- fibrosis1$time

# Quality control between the two arrays
par(mfrow=c(1,2))
fibrosis_diff <- (fibrosis1[, 3:ncol(fibrosis1)]-fibrosis2[, 3:ncol(fibrosis2)])
hist(as.matrix(abs(fibrosis_diff)), breaks=seq(from=0, to=5000, by=25),
     ylim=c(0,10), col="gray", 
     main="All probes",
     xlab="Difference fibrosis panels")
# there are some probes with very high difference. Prelimanary analyis shows this are the repeats
# of Mmp10. 
tmp <- fibrosis_diff
tmp[abs(tmp)<10] <- NA
tmp
# Removing Mmp10
hist(as.matrix(abs(fibrosis_diff[, colnames(fibrosis_diff)!="Mmp10"])), breaks=seq(from=0, to=5000, by=25),
     ylim=c(0,10), col="gray",
     main="All probes - Mmp10",
     xlab="Difference fibrosis panels")
par(mfrow=c(1,2))

# Remove the Mmp10 probe
fibrosis <- fibrosis[, colnames(fibrosis)!="Mmp10"]
rm(fibrosis1, fibrosis2, fibrosis_diff, tmp)

# subset of data is prepared and merged
d <- list()
d$gldh <- data.frame(histology[, c("sid_GLDH", "GLDH")])
names(d$gldh) <- c("sid", "GLDH")
d$alt <- data.frame(histology[, c("sid_ALT", "ALT")])
names(d$alt) <- c("sid", "ALT")
d$bilirubin <- data.frame(histology[, c("sid_bilirubin", "bilirubin")])
names(d$bilirubin) <- c("sid", "bilirubin")
d$albumin <- data.frame(histology[, c("sid_albumin", "albumin")])
names(d$albumin) <- c("sid", "albumin")
d$hc <- data.frame(histology[, c("sid_BrdU_HC", "BrdU_HC")])
names(d$hc) <- c("sid", "BrdU_HC")
d$nhc <- data.frame(histology[, c("sid_BrdU_NHC", "BrdU_NHC")])
names(d$nhc) <- c("sid", "BrdU_NHC")
d$kupffer <- data.frame(histology[, c("sid_BrdU_Kupffer", "BrdU_Kupffer")])
names(d$kupffer) <- c("sid", "BrdU_Kupffer")
d$bec <- data.frame(histology[, c("sid_BrdU_BEC", "BrdU_BEC")])
names(d$bec) <- c("sid", "BrdU_BEC")
d$siriusRed <- data.frame(histology[, c("sid_BrdU_SiriusRed", "BrdU_SiriusRed")])
names(d$siriusRed) <- c("sid", "BrdU_SirirusRed")
d$bileInfarcts <- data.frame(histology[, c("sid_bileInfarcts", "bileInfarcts")])
names(d$bileInfarcts) <- c("sid", "bileInfarcts")
summary(d)

# Merge the histological & antibody datasets on sample ids
tmp <- merge(d$gldh, d$alt, by="sid")
tmp <- merge(tmp, d$bilirubin, by="sid")
tmp <- merge(tmp, d$albumin, by="sid")
biochemistry <- tmp

tmp <- merge(d$hc, d$nhc, by="sid")
tmp <- merge(tmp, d$kupffer, by="sid")
tmp <- merge(tmp, d$bec, by="sid")
tmp <- merge(tmp, d$siriusRed, by="sid")
tmp <- merge(tmp, d$bileInfarcts, by="sid")
histology.processed <- tmp

# remove pid from antibodies
antibodies <- subset(antibodies, select = -c(pid) )

# Merge all data and create an overview of the factors
tmp <- merge(adme, cytokines, by = c('sid', 'time'))
tmp <- merge(tmp, fibrosis, by = c('sid', 'time'))
tmp <- merge(tmp, biochemistry, by=c('sid'))
tmp <- merge(tmp, histology.processed, by=c('sid'))
tmp <- merge(tmp, antibodies, by=c('sid'))
data <- tmp

# Type of factor
ftype <- c(rep("GE_ADME", ncol(adme)-2),
           rep("GE_Cytokines", ncol(cytokines)-2),
           rep("GE_Fibrosis", ncol(fibrosis)-2),
           rep("Biochemistry", ncol(biochemistry)-1),
           rep("Histology", ncol(histology.processed)-1),
           rep("Antibodies", ncol(antibodies)-1))
ftype.short <- c(rep("", ncol(adme)-2),
           rep("", ncol(cytokines)-2),
           rep("", ncol(fibrosis)-2),
           rep("B", ncol(biochemistry)-1),
           rep("H", ncol(histology.processed)-1),
           rep("A", ncol(antibodies)-1))

# Create the factor information
BDLfactors <- data.frame(id=colnames(data)[3:ncol(data)], ftype=ftype, ftype.short=ftype.short)


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
rm(tmp, d, adme, antibodies, cytokines, fibrosis, histology, histology.processed, biochemistry)


BDLdata <- data
BDLsamples <- samples
# save processed data
save(BDLdata, file="data/BDLdata.RData")
save(BDLsamples, file="data/BDLsamples.RData")
save(BDLfactors, file="data/BDLfactors.RData")
rm(data, samples, BDLdata, BDLsamples)


# information for probe mapping ------------------------------------------------
BDLprobes <- read.csv('inst/extdata/probe_mapping.csv', stringsAsFactors=FALSE, 
                      sep="\t")
save(BDLprobes, file="data/BDLprobes.RData")
