##############################################################################################
#
#    BDL - Analysis
#
#    Matthias Koenig
#    2015-08-03
#
#    This script contains the complete statistical analysis for the publication
#    Pathobiochemical signatures of cholestatic liver disease in bile duct ligated
#    mice (BMC Systems Biology).
#
#    Uses the BDLanalysis package to load all data and helper functions.
#
##############################################################################################

rm(list=ls())
require('BDLanalysis')
require('reshape')
require("calibrate")




help(BDLanalysis)

# Preprocess dataset
setwd("/home/mkoenig/git/bdl-analysis/BDLanalysis") 
source('data-raw/process.R')

# overview datasets
data()
# load the data sets 
data(BDLdata)
data(BDLsamples)
factors <- colnames(BDLdata)


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


# Mean of factors for the individual time points
tmp <- bdl_mean_data(data)
dmean <- tmp$dmean
dmean.time <- tmp$dmean.time
rm(tmp)

head(dmean)
head(dmean.time)

# Data list
data_list <- bdl_data_matrices(data, time_pts=levels(time))



# plot options
# @export
# options = list(width=1600, height=800, res=150)



