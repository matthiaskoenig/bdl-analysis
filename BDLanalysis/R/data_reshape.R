# -----------------------------------------------------------------------------------
# Data reshaping of raw data. 
# This consists for instance in the caluclation of mean data,
# or the reshaping into data matrices.
# -----------------------------------------------------------------------------------

#' Calculate mean data average over time points.
#' 
#' Use the time point repeats for averaging. The time points are defined
#' in the sample definitions
#' @export
bdl_mean_data <- function(data, samples){  
  data2 <- data
  data2$time <- samples$time
  # rm NA in mean calculation
  dmean <- aggregate(data2, list(data2$time), FUN=mean, na.rm=TRUE)
  time <- dmean$time
  dmean <- subset(dmean, select = -c(time, Group.1))
  rownames(dmean) <- levels(samples$time_fac)
  return(dmean)
}


#' Transform BDL data frame into list of matrices.
#' 
#' @export
bdl_matrix_data <- function(data, samples, Nrepeats=5){
  time_pts <- levels(samples$time_fac)
  data_list <- list()
  for (name in names(data)){
    # important to fill in the right order !
    data_list[[name]] <- matrix(data[[name]], nrow=length(time_pts), ncol=Nrepeats, byrow=TRUE)
    colnames(data_list[[name]]) <- paste('R', 1:Nrepeats, sep="")
    rownames(data_list[[name]]) <- time_pts
  }
  names(data_list) <- names(data)
  return(data_list)
}
