#---------------------------------------------------------------------------------------------
# Read & preprocess data
#
# The individual raw data files are loaded and combined in a single data frame.
# Factor and timecourse information is added.
# All processed datasets are stored as Rdata files in combination with the packages.
#
# This should only be run to recreate the processed data files after the raw data
# files changed.
# 
# TODO: preprocess all the data into Rdata files which can be loaded in the analysis scripts.
#---------------------------------------------------------------------------------------------

samples <- read.csv('data-raw/samples.csv', sep="\t")
histology <- read.csv('data-raw/histology.csv', sep="\t")
adme <- read.csv('data-raw/fluidigm_ADME.csv', sep="\t")
cytokines <- read.csv('data-raw/fluidigm_cytokines.csv', sep="\t")
fibrosis1 <- read.csv('data-raw/fluidigm_fibrosis_01.csv', sep="\t")
fibrosis2 <- read.csv('data-raw/fluidigm_fibrosis_02.csv', sep="\t")
antibodies <- read.csv('data-raw/antibodies.csv', sep="\t")

# Two repeats of the fibrosis fluidigm chip were measured. 
# The mean value of the two chips is used for analysis.
fibrosis <- (fibrosis1 + fibrosis2) / 2
fibrosis$time <- fibrosis1$time
rm(fibrosis1, fibrosis2)

# subset of data is prepared and merged
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


BDLdata <- data
BDLsamples <- samples
# save processed data
save(BDLdata, file = "data/BDLdata.rdata")
save(BDLsamples, file = "data/BDLsamples.rdata")
rm(data, samples, BDLdata, BDLsamples)
