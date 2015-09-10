#    BDL - Regression Analysis
This repository contains all data, information and the complete statistical analysis for the publication
**Pathobiochemical signatures of cholestatic liver disease in bile duct ligated mice** (BMC Systems Biology).
The analysis was performed in `R` and is available via the package `BDLanalysis` and the package vignette `Koenig_BDL_analysis.Rmd`.

Pathobiochemical signatures of cholestatic liver disease in bile duct ligated (BDL) mice were analysed in a comprehensive data set of serum markers, (immuno-)histological parameters and transcript profiles at 8 time points after BDL in mice, comprising different stages of the disease (> 6000 data points).

The main steps of the analysis comprise

* Explorative data anaysis
* Dimension reduction via ANOVA
* Correlation analysis
* Hierarchical clustering
* Decision trees

**author**: Matthias Koenig  
**version**: 0.3  

## Run analysis
For running the analysis checkout the git repository and open the `BDLanalysis.Rproj` with `RStudio`. You can than build the `BDLanalysis` package and recreate the package vignette `Koenig_BDL_analysis.Rmd` with `knitr`.
