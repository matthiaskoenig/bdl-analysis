#' BDLanalysis
#'
#' @name BDLanalysis
#' @docType package
NULL


#' Time course data from BDL mice.
#'
#' Complete data set of all time course data for all factors. The sample
#' description is available in BDLsamples.
#'
#' @format A data frame with with 40 rows and 154 variables: \code{Ppara}, \code{Cyp3a11},
#'   \code{Cyp24a1}, ...).
#' @name BDLdata
#' @docType data   
NULL


#' Sample description corresponding to BDLdata.
#'
#' Sample description for the BDLdata factor data set. Here the time point and
#' repeat information is defined.
#'
#' @format A data frame with with 40 rows and 9 variable providing additional information for the samples: \code{time}, \code{time_fac}, \code{time_point}, \code{repeats} ...).
#' @name BDLsamples
#' @docType data
NULL

#' Fluidigm probe mapping to UniProt.
#'
#' Maps the gene identifier to UniProt identifiers for additional information
#' retrieval. This information is used for customizing the plots.
#'
#' @format A data frame with with 40 rows and 9 variable providing additional information for the samples: \code{time}, \code{time_fac}, \code{time_point}, \code{repeats} ...).
#' @name BDLprobes
#' @docType data
NULL