##############################################################################################
#    Data reshaping of raw data
##############################################################################################

#' Calculate mean data average over time points.
#' 
#' Use the time point repeats for averaging. The time points are defined
#' in the sample definitions
#' @export
bdl_mean_data <- function(data, samples){  
  dmean <- bdl_fun_data(data=data, samples=samples, fun=mean)
  return(dmean)
}


#' Calculate sd data over time points.
#' 
#' Use the time point repeats for calculation of sd. The time points are defined
#' in the sample definitions
#' @export
bdl_sd_data <- function(data, samples){  
  dsd <- bdl_fun_data(data=data, samples=samples, fun=sd)
  return(dsd)
}


#' Calculate function fun over time points.
#' 
#' The function is used to calculate the mean and sd on the BDL data set.
#' @export
bdl_fun_data <- function(data, samples, fun){  
  data2 <- data
  data2$time <- samples$time
  # rm NA in mean calculation
  dfun <- aggregate(data2, list(data2$time), FUN=fun, na.rm=TRUE)
  time <- dfun$time
  dfun <- subset(dfun, select = -c(time, Group.1))
  rownames(dfun) <- levels(samples$time_fac)
  return(dfun)
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
