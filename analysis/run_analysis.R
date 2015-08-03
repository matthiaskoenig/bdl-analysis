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
setwd("/home/mkoenig/git/bdl-analysis/analysis") # set this to your analysis folder

# plot options
# @export
# options = list(width=1600, height=800, res=150)



# data preprocessing
source('data_processing.R')