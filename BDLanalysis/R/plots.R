##############################################################################################
#    Helper functions to generate plots for certain data sets.
##############################################################################################


#' Create a single factor plot.
#'
#' Plot the time course of a single factor, i.e. all individual data points.
#' @export
plot_single_factor <- function(name){
  dA <- BDLdata[, name]
  
  # search if probe info is available
  info <- ProbeInformation(geneId=name)  
  
  par(mfrow=c(1,2))
  
  # [A] plot against time
  plot(BDLsamples$time_fac, dA, at=sort(as.numeric(levels(as.factor(BDLsamples$time)))), col="blue",
       xlab="time [h]", ylab=name, main=name,
       ylim=c(0, max(dA, na.rm=TRUE)*1.1 ))
  
  points(BDLsamples$time, dA, col="black")
  points(BDLsamples$time, dA, col=rgb(0,0,1,0.6), pch=16)
  if (!is.null(info$Protein.name)){
    text(x=140, y=max(dA, na.rm=TRUE)*1.08, 
         labels=info$Protein.name, cex=0.8)
  }
  # add mean
  mean.time <- as.numeric(levels(as.factor(BDLsamples$time)))
  Nt <- length(mean.time)
  points(mean.time, BDLmean[, name], col="red", pch=15)
  lines(mean.time, BDLmean[, name], col="red")
  
  # [B] plot as factor (non-equidistant time points)
  plot(BDLsamples$time_fac, dA, xlab="time", ylab=name, main=name, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1))
  points(BDLsamples$time_fac, dA, col="black")
  points(BDLsamples$time_fac, dA, col=rgb(0,0,1,0.6), pch=16)
  # add mean
  points(1:Nt, BDLmean[, name], col="red", pch=15)
  lines(1:Nt, BDLmean[, name], col="red")
  
  par(mfrow=c(1,1))
}

#' Creates single factor plots of all factors.
#' 
#' Iterates over all factor variables in the BDLdata frame and creates
#' the individual factor plots.
#' @export 
plot_all_factors <- function(){
  
  factors <- colnames(BDLdata)
  Nf = length(factors)
  for (k in 1:Nf){
    cat(sprintf("%s / %s\n", k, Nf))
    name <- factors[k]
    # fname <- paste("../results/factors/", sprintf("%03d", k), "_", name, ".png", sep="")
    # png(filename=fname, width=options$width, height=options$height, res=options$res)
    plot_single_factor(name)
    # dev.off()
  }  
}

#' Plot of a single factor.
#' 
#' Plots the mean time course and individual samples.
#' @export
f_single_plot <- function(name_A){
  dA <- data[[name_A]]
  dA.mean <- dmean[[name_A]]
  plot(samples$time_fac, dA, xlab="time", ylab=name_A, main=name_A, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1))
  points(samples$time_fac, dA, col="black")
  points(samples$time_fac, dA, col=rgb(0,0,1,0.6), pch=16)
  # plot the repeat number
  dA.text <- dA
  dA[is.na(dA)] <- -1
  # library(calibrate)
  textxy(samples$time_point, dA, samples$sid, col="black", cex=0.6)
  
  points(1:nrow(dmean), dA.mean, col="red", pch=15)
  lines(1:nrow(dmean), dA.mean, col="red")
}

#' Correlation plot between two factors
#' 
#' @export
f_cor_pair_plot <- function(name_A, name_B, single_plots=TRUE){
  if (single_plots){
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
  }
  # correlation plot
  dA <- data[[name_A]]
  dB <- data[[name_B]]
  # dA <- dA[!is.na(dA)]
  # d <- dA[!is.na(dA)]
  dA.mean <- dmean[[name_A]]
  dB.mean <- dmean[[name_B]]
  
  value.max <- 1.1* c(max(dA), max(dB))
  plot(dA, dB, xlim=c(0, value.max[1]), ylim=c(0, value.max[2]), 
       main=sprintf("%s ~ %s", name_A, name_B),
       xlab=name_A, ylab=name_B)
  textxy(dA, dB, samples$time)
  points(dA.mean, dB.mean, pch=16, col=rgb(0,0,1,0.8) )
  lines(dA.mean, dB.mean, col=rgb(0,0,1,0.8) )
  textxy(dA.mean, dB.mean, dmean.time, col="blue")
  abline(a=0, b=mean(dB.mean)/mean(dA.mean), col="darkgray")
  
  text(x=0.05*value.max[1], y=0.9*value.max[2], 
       labels=sprintf("S: %1.3f\nS mean: %1.3f\nP: %1.3f\nP mean: %1.3f", 
                      cor(dA, dB, method="spearman", use="pairwise.complete.obs"),
                      cor(dA.mean, dB.mean, method="spearman", use="pairwise.complete.obs"),
                      cor(dA, dB, method="pearson", use="pairwise.complete.obs"),
                      cor(dA.mean, dB.mean, method="pearson", use="pairwise.complete.obs")))
  
  # single plots
  if (single_plots){
    f_single_plot(name_A)
    f_single_plot(name_B)
    layout(matrix(c(1), 1, 1, byrow = TRUE))
  }
}


#' Fluidigm probe annotation
#'
#' Helper function with information for fluidigm probes.
#' @export
ProbeInformation <- function(geneId){
  info <- list()
  idx <- which(BDLprobes$Gene==geneId)
  if (!is.null(idx) & length(idx)>0){
    info$Protein.name <- BDLprobes$Protein.names[idx[1]]
    info$Chip <- BDLprobes$Chip[idx[1]]
    info$Entry <- BDLprobes$Entry[idx[1]]
    info$Gene.name <- BDLprobes$Gene.names[idx[1]]
  } else {
    info$Protein.name <- NULL
  }
  return(info)
}
