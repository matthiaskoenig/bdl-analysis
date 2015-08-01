##############################################################################################
#
#    BDL - Analysis
#
#    Matthias Koenig
#    2015-08-01
#
#    This script contains the complete statistical analysis for the publication
#    Pathobiochemical signatures of cholestatic liver disease in bile duct ligated
#    mice (BMC Systems Biology).
#
#    The necessary requirements can be installed with install.R.
#
##############################################################################################

rm(list=ls())
setwd("/home/mkoenig/git/bdl-analysis/analysis") # set this to your analysis folder

# data preprocessing
source('data_processing.R')