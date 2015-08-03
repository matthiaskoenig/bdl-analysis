# -----------------------------------------------------------------------------------
# Data reshaping of raw data. 
# This consists for instance in the caluclation of mean data,
# or the reshaping into data matrices.
# -----------------------------------------------------------------------------------

#' Calculation of mean data based on provided time factor.
#' 
#' @export
bdl_mean_data <- function(data){  
  data2 <- data
  data2$time <- samples$time
  # rm NA in mean calculation
  dmean <- aggregate(data2, list(data2$time), FUN=mean, na.rm=TRUE)
  dmean.time <- dmean$time
  dmean <- subset(dmean, select = -c(time, Group.1))
  return( list(dmean=dmean, dmean.time=dmean.time) )
}


#' List of data matrices.
#' 
#' @export
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
