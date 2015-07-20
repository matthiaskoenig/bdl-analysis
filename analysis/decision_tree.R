##############################################################################################
#
#    BDL - Decision Tree
#
#    Matthias Koenig
#    2015-07-20
#
##############################################################################################

rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis")

# Load data for prediction
file.data <- file.path("..", "data", "bdl-data.Rdata")
load(file=file.data)

# The class to predict is the actual time point
samples$time_fac

# The possible predictors are all factors in the data set
head(data)