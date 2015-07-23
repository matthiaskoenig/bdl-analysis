##############################################################################################
#
#    BDL - Regression Analysis
#
#    Matthias Koenig
#    2015-07-21
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
# install.packages("calibrate")
require("calibrate")


#---------------------------------------------
# Read & preprocess data
#---------------------------------------------
# TODO: better merge of data
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

#---------------------------------------------
# Dimension reduction via ANOVA
#---------------------------------------------
# A one-way analyis of variance (ANOVA) was applied to filter genes showing
# significant (padj<0.05) up- or down-regulation during the time course, using the Bonferonni
# step-down procedure to correct for any artificial p-value inflation.

# http://www.r-tutor.com/elementary-statistics/analysis-variance/completely-randomized-design

# perform anova for the factor
single_factor_anova <- function(name){
  # data matrix
  mat1 <- t(data_list[[name]])
  colnames(mat1) <- levels(samples$time_fac)
  
  # Concatenate the data rows of df1 into a single vector r .
  r = c(t(as.matrix(mat1))) # response data 
  
  # Assign new variables for the treatment levels and number of observations.
  f = levels(samples$time_fac)   # treatment levels 
  k = 8                          # number of treatment levels 
  n = 5                          # observations per treatment 
  
  # Create a vector of treatment factors that corresponds to each element of r in step 3 with the gl function.
  tm = gl(k, 1, n*k, factor(f))   # matching treatments 
  
  # Apply the function aov to a formula that describes the response r by the treatment factor tm.
  # Fit an analysis of variance model
  av = aov(r ~ tm) 
  
  # Print out the ANOVA table with the summary function. 
  # summary(av)
  return(av)
}

# Perform sanova for all factors.
# The unadjusted p-values are returned.
all_factor_anova <- function(){
  df.anova <- data.frame(factors)
  df.anova$p.value <- NA
  
  for (k in 1:nrow(df.anova)){
    name <- factors[[k]]
    av <- single_factor_anova(name)
    p.value <- summary(av)[[1]][["Pr(>F)"]][[1]]
    df.anova[k, "p.value"] <- p.value
  }
  return(df.anova)
}

# creates significant codes for given p.value.
significant_code <- function(p.value){
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  sig = " "
  if (p.value <= 0.001){
    sig = "***"
  } else if (p.value <= 0.01){
    sig = "**"
  } else if (p.value <=0.05){
    sig = "*"
  } else if (p.value <=0.1){
    sig = "."
  } else if (p.value <=1){
      sig = " "
  }
  return (sig)
}
# calculate anova
df.anova <- all_factor_anova()
df.anova$sig <- sapply(df.anova$p.value, significant_code)

# Adjust the p-values for multiple testing
# Given a set of p-values, returns p-values adjusted using one of several methods.
# The Bonferroni, Holm, Hochberg, Hommel are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions. 
# Using Holm correction with number of tests.
# Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65–70.
df.anova$p.holm <- p.adjust(df.anova$p.value, method ="holm" , n = length(df.anova$p.value))
df.anova$sig.holm <- sapply(df.anova$p.holm, significant_code)

df.anova

# save the ordered results of the ANOVA test
df.anova.ordered <- df.anova[with(df.anova, order(p.holm)), ]
summary(df.anova.ordered)
write.table(df.anova.ordered, file="../results/factor_anova.csv", sep="\t", quote=FALSE)

file.anova <- file.path("..", "data", "bdl-anova.Rdata")
save(df.anova, file=file.anova)

# if filtering is selected the subset of factors being different 
# within the timecourse are selected
filter_by_anova = TRUE
if (filter_by_anova){
  # acceptance level
  p.accept = 0.05
  
  # how many rejected by adjusted p-value
  table(df.anova$p.holm>=p.accept) # 64 rejected / 90 accepted
  table(df.anova$p.value>=p.accept) # 19 rejected / 135 accepted
  
  # the accepted subset
  f.accept = df.anova$p.holm<p.accept
  
  # overwrite the full data set
  data.full <- data
  data <- data[, f.accept]
  factors.full <- factors
  factors <- factors[f.accept]
  dmean.full <- dmean
  dmean <- dmean[, f.accept]  
}

# plot of the data subset which is used for the correlation analysis
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")) 
heatmap.2(t(as.matrix(data)), col=col2(100), scale="row", Rowv=NULL, Colv=NULL,
          key=TRUE, trace="none", cexRow=0.5, keysize=0.8)

#---------------------------------------------
# Correlation analysis
#---------------------------------------------
# Calculation of standard correlation matrices based on Spearman
# and Pearson correlaton coefficients. 
# Hierarchical clustering is performed using complete linkage.
# The cor.pearson and cor.spearman are used as comparison models.

# Do the analysis either on the full dataset or the filtered subset




# correlation matrix
cor.pearson <- cor(data, method="pearson", use="pairwise.complete.obs")
cor.spearman <- cor(data, method="spearman", use="pairwise.complete.obs")
cor.pearson.fil <- cor(data.fil, method="pearson", use="pairwise.complete.obs")
cor.spearman.fil <- cor(data.fil, method="spearman", use="pairwise.complete.obs")


library(corrplot)
options <- list(width=2000, height=2000, res=200)

# helper function for correlation plot
f_corrplot <- function(name, data, order){
  fname <- sprintf("%s_%s.png", name, order)
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  png(filename=sprintf("../results/%s", fname), width=options$width, height=options$height, res=options$res)
  corrplot(data, order=order, hclust.method="complete", method="color", type="full", 
           tl.cex=0.3, tl.col="black", col=col2(10))
  dev.off()
}

# Spearman & Pearson correlation matrices
f_corrplot("cor.spearman", data=cor.spearman, order="original")
f_corrplot("cor.spearman", data=cor.spearman, order="hclust")
f_corrplot("cor.pearson", data=cor.pearson, order="original")
f_corrplot("cor.pearson", data=cor.pearson, order="hclust")
f_corrplot("cor.spearman.fil", data=cor.spearman.fil, order="original")
f_corrplot("cor.spearman.fil", data=cor.spearman.fil, order="hclust")
f_corrplot("cor.pearson.fil", data=cor.pearson.fil, order="original")
f_corrplot("cor.pearson.fil", data=cor.pearson.fil, order="hclust")



# --- YR1, YS1, YR2, YS2 correlation ----------
# Calculation of time-course based correlation measurements, namely ys1, ys2, yr1, yr2
# and the respective adaptions to the underlying datasets, i.e. using multiple
# repeat data with every time point measurment coming from an individual
# sample.
# In ys1, ys2, yr1, yr2 all calculations are performed on the mean time course data.
# The classical correlation components are replaced with the correlations calculated
# on the individual data points.

source("ys1_yr1.R")        # definition of ys1, ys2, yr1, yr2
source("ys1_yr1_tools.R")  # definition of ys1 and yr1

# calculation of ys1 and yr1 for mean data
w <- list(w1=0.5, w2=0.25, w3=0.25)
res.ys1 <- ys1.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1 <- res.ys1$value
cor.ys1.scaled <- 2*(cor.ys1-0.5)
res.ys1.fil <- ys1.df(dmean.fil, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys1.fil <- res.ys1.fil$value
cor.ys1.fil.scaled <- 2*(cor.ys1.fil-0.5)

res.ys2 <- ys2.df(dmean, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2 <- res.ys2$value
cor.ys2.scaled <- 2*(cor.ys2-0.5)
res.ys2.fil <- ys2.df(dmean.fil, dmean.time, w1=w$w1, w2=w$w2, w3=w$w3, use="pairwise.complete.obs")
cor.ys2.fil <- res.ys2.fil$value
cor.ys2.fil.scaled <- 2*(cor.ys2.fil-0.5)


f_corrplot("cor.ys1", data=cor.ys1, order="original")
f_corrplot("cor.ys1", data=cor.ys1, order="hclust")
f_corrplot("cor.ys1.fil", data=cor.ys1.fil, order="original")
f_corrplot("cor.ys1.fil", data=cor.ys1.fil, order="hclust")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="original")
f_corrplot("cor.ys1.A", data=res.ys1$A, order="hclust")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="original")
f_corrplot("cor.ys1.M", data=res.ys1$M, order="hclust")

f_corrplot("cor.ys2", data=cor.ys2, order="original")
f_corrplot("cor.ys2", data=cor.ys2, order="hclust")
f_corrplot("cor.ys2.fil", data=cor.ys1.fil, order="original")
f_corrplot("cor.ys2.fil", data=cor.ys1.fil, order="hclust")

f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="original")
f_corrplot("cor.ys2.A_star", data=res.ys2$A_star, order="hclust")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="original")
f_corrplot("cor.ys2.M_star", data=res.ys2$M_star, order="hclust")

f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="original")
f_corrplot("cor.ys1.scaled", data=cor.ys1.scaled, order="hclust")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="original")
f_corrplot("cor.ys2.scaled", data=cor.ys2.scaled, order="hclust")
f_corrplot("cor.ys1.fil.scaled", data=cor.ys1.fil.scaled, order="original")
f_corrplot("cor.ys1.fil.scaled", data=cor.ys1.fil.scaled, order="hclust")
f_corrplot("cor.ys2.fil.scaled", data=cor.ys2.fil.scaled, order="original")
f_corrplot("cor.ys2.fil.scaled", data=cor.ys2.fil.scaled, order="hclust")

# Pearson & spearman correlation on full dataset as replacement for the 
# mean Pearson/Spearman in ys1, ys2, yr1, yr2
# Now calculate the scores with full correlation
cor.S_star <- ( cor(data, method="spearman", use="pairwise.complete.obs") + 1 )/2
cor.R_star <- ( cor(data, method="pearson", use="pairwise.complete.obs") + 1 )/2
cor.S_star.fil <- ( cor(data.fil, method="spearman", use="pairwise.complete.obs") + 1 )/2
cor.R_star.fil <- ( cor(data.fil, method="pearson", use="pairwise.complete.obs") + 1 )/2


w <- list(w1=0.5, w2=0.25, w3=0.25)
cor.ys1_full <- w$w1*cor.S_star + w$w2*res.ys1$A + w$w3*res.ys1$M
cor.ys1_full.scaled <- 2*(cor.ys1_full-0.5)
cor.ys2_full <- w$w1*cor.S_star + w$w2*res.ys2$A_star + w$w3*res.ys2$M_star
cor.ys2_full.scaled <- 2*(cor.ys2_full-0.5)
# considering slope and time difference
cor.ys3_full <- w$w1*cor.S_star + w$w2*res.ys2$A_star2 + w$w3*res.ys2$M_star2
cor.ys3_full.scaled <- 2*(cor.ys3_full-0.5)

cor.yr1_full <- w$w1*cor.R_star + w$w2*res.ys1$A + w$w3*res.ys1$M
cor.yr1_full.scaled <- 2*(cor.yr1_full-0.5)
cor.yr2_full <- w$w1*cor.R_star + w$w2*res.ys2$A_star + w$w3*res.ys2$M_star
cor.yr2_full.scaled <- 2*(cor.yr2_full-0.5)
# considering slope and time difference
cor.yr3_full <- w$w1*cor.R_star + w$w2*res.ys2$A_star2 + w$w3*res.ys2$M_star2
cor.yr3_full.scaled <- 2*(cor.yr3_full-0.5)

f_corrplot("cor.ys1_full.scaled", data=cor.ys1_full.scaled, order="original")
f_corrplot("cor.ys1_full.scaled", data=cor.ys1_full.scaled, order="hclust")
f_corrplot("cor.ys2_full.scaled", data=cor.ys2_full.scaled, order="original")
f_corrplot("cor.ys2_full.scaled", data=cor.ys2_full.scaled, order="hclust")
f_corrplot("cor.ys3_full.scaled", data=cor.ys3_full.scaled, order="original")
f_corrplot("cor.ys3_full.scaled", data=cor.ys3_full.scaled, order="hclust")

f_corrplot("cor.yr1_full.scaled", data=cor.yr1_full.scaled, order="original")
f_corrplot("cor.yr1_full.scaled", data=cor.yr1_full.scaled, order="hclust")
f_corrplot("cor.yr2_full.scaled", data=cor.yr2_full.scaled, order="original")
f_corrplot("cor.yr2_full.scaled", data=cor.yr2_full.scaled, order="hclust")
f_corrplot("cor.yr3_full.scaled", data=cor.yr3_full.scaled, order="original")
f_corrplot("cor.yr3_full.scaled", data=cor.yr3_full.scaled, order="hclust")

# on filtered data
cor.ys1_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys1.fil$A + w$w3*res.ys1.fil$M
cor.ys1_full.fil.scaled <- 2*(cor.ys1_full.fil-0.5)
cor.ys2_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys2.fil$A_star + w$w3*res.ys2.fil$M_star
cor.ys2_full.fil.scaled <- 2*(cor.ys2_full.fil-0.5)
# considering slope and time difference
cor.ys3.fil <- w$w1*res.ys2.fil$S_star + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.ys3.fil.scaled <- 2*(cor.ys3.fil-0.5)
cor.ys3_full.fil <- w$w1*cor.S_star.fil + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.ys3_full.fil.scaled <- 2*(cor.ys3_full.fil-0.5)

cor.yr1_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys1.fil$A + w$w3*res.ys1.fil$M
cor.yr1_full.fil.scaled <- 2*(cor.yr1_full.fil-0.5)

cor.yr2_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys2.fil$A_star + w$w3*res.ys2.fil$M_star
cor.yr2_full.fil.scaled <- 2*(cor.yr2_full.fil-0.5)
# considering slope and time difference
cor.yr3_full.fil <- w$w1*cor.R_star.fil + w$w2*res.ys2.fil$A_star2 + w$w3*res.ys2.fil$M_star2
cor.yr3_full.fil.scaled <- 2*(cor.yr3_full.fil-0.5)

f_corrplot("cor.ys1_full.fil.scaled", data=cor.ys1_full.fil.scaled, order="hclust")
f_corrplot("cor.ys2_full.fil.scaled", data=cor.ys2_full.fil.scaled, order="hclust")
f_corrplot("cor.ys3_full.fil.scaled", data=cor.ys3_full.fil.scaled, order="hclust")

f_corrplot("cor.yr1_full.fil.scaled", data=cor.yr1_full.fil.scaled, order="hclust")
f_corrplot("cor.yr2_full.fil.scaled", data=cor.yr2_full.fil.scaled, order="hclust")
f_corrplot("cor.yr3_full.fil.scaled", data=cor.yr3_full.fil.scaled, order="hclust")


# name_A <- "Actb"
# name_B <- "Actb.x"
# extreme example where the difference between pearson and spearman matters
name_A <- "Nos2"
name_B <- "Cxcl15"
f_cor_pair_plot("Nos2", "Cxcl15")

f_cor_pair_plot("albumin", "Cyp2b10")

#---------------------------------------------
# Hierarchical clustering
#---------------------------------------------
# factor 1, factor 2, pearson, spearman, ys1_mean, ys1, ...
# -> plot the clusters and cluster members

# Perform hierarchical clustering on matrix
# The hclust function in R uses the complete linkage method for hierarchical clustering by default.
# This particular clustering method defines the cluster distance between two clusters to be the 
# maximum distance between their individual components.

# apply hirarchical clustering 
method = "yr3"
if (identical(method, "ys1")){
  cor.cluster <- cor.ys1_full.scaled  
}else if (identical(method, "ys2")){
  cor.cluster <- cor.ys2_full.scaled  
}else if (identical(method, "ys3")){
  cor.cluster <- cor.ys3_full.scaled  
}else if (identical(method, "yr1")){
  cor.cluster <- cor.yr1_full.scaled  
}else if (identical(method, "yr2")){
  cor.cluster <- cor.yr2_full.scaled  
}else if (identical(method, "yr3")){
  cor.cluster <- cor.yr3_full.scaled  
}

method = "yr1"
if (identical(method, "ys1")){
  cor.cluster <- cor.ys1_full.fil.scaled  
}else if (identical(method, "ys2")){
  cor.cluster <- cor.ys2_full.fil.scaled  
}else if (identical(method, "ys3")){
  cor.cluster <- cor.ys3_full.fil.scaled  
}else if (identical(method, "yr1")){
  cor.cluster <- cor.yr1_full.fil.scaled  
}else if (identical(method, "yr2")){
  cor.cluster <- cor.yr2_full.fil.scaled  
}else if (identical(method, "yr3")){
  cor.cluster <- cor.yr3_full.fil.scaled  
}


test <- dist(cor.cluster)

hc <- hclust(dist(cor.cluster)) 
corrplot(cor.cluster[hc$order, hc$order], order="original", method="square", type="full",tl.cex=0.3, tl.col="black")

# plot(hc)               # plot the dendrogram 
plot(hc, hang=-1)               # plot the dendrogram 

# cut tree into clusters
Ngroups = 8
rect.hclust(hc, k=Ngroups)
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)
groups.hc.order <- groups[hc$order]

dev.off()

# plot the groups
options$height = 3000
options$width = 3000

# Plot the clusters
for (k in 1:Ngroups){
  fname <- sprintf("%s_cluster_%s.png", method, k)
  png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  par(mfrow=c(N,N))
  for (name in names(g)){
    f_single_plot(name_A=name) 
  }
  par(mfrow=c(1,1))  
  dev.off()
}

#install.packages('matrixStats')
library('matrixStats')

# mean plots for clusters 
# TODO: variance to the cluster
f_normalize_centering <- function(a){
  a.norm <- (a - mean(a))/(max(a, na.rm=TRUE) - min(a, na.rm=TRUE))
  return(a.norm)
}
options$height = 1600
options$width = 1600
fname <- sprintf("%s_cluster_overview.png", method)
png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
par(mfrow=c(ceiling(sqrt(Ngroups)),ceiling(sqrt(Ngroups))))
steps = 1:8
for (k in 1:Ngroups){
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  dgroup <- dmean[names(g)]

  # centralize and normalize columns, i.e. the individual factors for comparison
  dgroup.norm <- apply(dgroup, 2, f_normalize_centering)
  
  # mean and sd for timepoints 
  g.mean <- rowMeans(dgroup.norm)
  g.sd <- rowSds(dgroup.norm)   # apply(dgroup.norm, 2, sd)
  
  # plot sd range
  plot(1, type="n", xlab="", ylab="", xlim=c(1, 8), ylim=c(-1, 1), main=sprintf("%s : Cluster %s", method, k))
  polygon(c(steps, rev(steps)), c(g.mean+g.sd, rev(g.mean-g.sd)),
          col = rgb(0.5,0.5,0.5,0.5), border = NA)
  
  # individual data
  for (name in names(g)){
    points(steps, dgroup.norm[, name], pch=16, col="black")
    lines(steps, dgroup.norm[, name], col=rgb(0.5,0.5,0.5, 0.8), lwd=1)
  }
  # mean over factors in cluster
  lines(steps, g.mean, col="blue", lwd=2)
}
par(mfrow=c(1,1))
dev.off()


# better dendrogram
dend1 <- as.dendrogram(hc)
# Get the package:
# install.packages("dendextend")
library(dendextend)

# Get the package:
cutree(dend1,h=70) # it now works on a dendrogram
# It is like using:
dendextend:::cutree.dendrogram(dend1,h=70)
# dend1 <- color_branches(dend1, k = 7)
dend1 <- color_labels(dend1, k=Ngroups)
plot(dend1)



# TODO: create the list of all pairwise correlations
# TODO: GO annotations of clusters


f_cor_pair_plot("Cebpa", "Cyp24a1")
  