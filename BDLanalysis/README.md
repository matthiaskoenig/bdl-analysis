#    BDL - Regression Analysis

author: Matthias Koenig  
date: 2015-08-03  
version: 0.1

This repository contains all data, information and the complete statistical analysis for the publication
**Pathobiochemical signatures of cholestatic liver disease in bile duct ligated mice** (BMC Systems Biology).
The analysis was performed in `R` and is available as packacke `BDLanalysis`.

Pathobiochemical signatures of cholestatic liver disease in bile duct ligated (BDL) mice were analysed 
in a comprehensive data set of serum markers, histological parameters and transcript profiles at 8 time points after bile duct ligation in mice, comprising different stages of the disease (> 6000 data points).

Analysis consists among others of
* filtering of relevant factors via ANOVA
* correlation analysis using correlation measure for time course data
* clustering of factors with hiearachical clustering methods
* time class prediction with decision trees

## Run analysis
For running the analysis the `BDLanalysis` package has to be installed.

Install it from github with:
```
devtools::install_github("matthiaskoenig/BDLanalysis")
```

