% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ys1_yr1.R
\name{ysr.matrices}
\alias{ysr.matrices}
\title{Calculate ys and yr component matrices given data frame.}
\usage{
ysr.matrices(data, time_pts, use = "pairwise.complete.obs")
}
\description{
Every column of the data.frame is a single measurement (or mean) of a single factor for
the provided time points.
Use this function to calculate the correlation matrix on the mean data.
}

