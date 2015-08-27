##############################################################################################
#    Helper functions to generate plots for certain data sets.
##############################################################################################


#' Create a single factor plot.
#'
#' Plot the time course of a single factor, i.e. all individual data points.
#' @export
plot_single_factor <- function(name, path=NULL, k=NULL){
  # get factor for name
  dA <- BDLdata[, name]
  # probe info for gene
  info <- ProbeInformation(geneId=name)  
  
  if (! is.null(path)){
    fname <- file.path(path, paste(sprintf("%03d", k), "_", name, ".png", sep=""))
    png(filename=fname, width=1600, height=800, res=150)
  }
  par(mfrow=c(1,2))
  title <- sprintf("%s (%s)", name, BDLfactors$ftype[BDLfactors$id==name])
  
  # [A] plot against time
  plot(BDLsamples$time_fac, dA, at=sort(as.numeric(levels(as.factor(BDLsamples$time)))), col="blue",
       xlab="time [h]", ylab=name, main=title, font.lab=2,
       xlim=c(-0.5, max(BDLsamples$time)*1.1),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1 ),
       cex.axis=0.8)
  
  points(BDLsamples$time, dA, col="black")
  points(BDLsamples$time, dA, col=rgb(0,0,1,0.6), pch=16)
  # labels
  require(calibrate)
  na.idx <- is.na(dA)
  textxy(BDLsamples$time[!na.idx], dA[!na.idx], BDLsamples$sid, col="black", cex=0.5)
  # protein name
  if (!is.null(info$Protein.name)){
    text(x=140, y=max(dA, na.rm=TRUE)*1.08, 
         labels=info$Protein.name, cex=0.8)
  }
  # mean
  Nt <- length(BDLmean.time)
  points(BDLmean.time, BDLmean[, name], col="red", pch=15)
  lines(BDLmean.time, BDLmean[, name], col="red")
  
  # [B] plot as factor (non-equidistant time points)
  plot(BDLsamples$time_fac, dA, xlab="time [class]", ylab=name, main=title, col=rgb(0.5,0.5,0.5, 0.4),
       font.lab=2, ylim=c(0, max(dA, na.rm=TRUE)*1.1), cex.axis=0.8)
  points(BDLsamples$time_fac, dA, col="black")
  points(BDLsamples$time_fac, dA, col=rgb(0,0,1,0.6), pch=16)
  # labels
  textxy(as.numeric(BDLsamples$time_fac[!na.idx]), dA[!na.idx], BDLsamples$sid, col="black", cex=0.5)
  # add mean
  points(1:Nt, BDLmean[, name], col="red", pch=15)
  lines(1:Nt, BDLmean[, name], col="red")
  
  legend("topleft", col=c("red", rgb(0,0,1,0.6)), 
         legend=c("time point mean", "single sample"), bty="n", 
         pch=c(15,16), cex=0.7)
  
  par(mfrow=c(1,1))
  if (!is.null(path)){
    dev.off()
  }
}

#' Creates single factor plots of all factors.
#' 
#' Iterates over all factor variables in the BDLdata frame and creates
#' the individual factor plots. BDLdata has to be available in the 
#' environment
#' @export 
plot_all_factors <- function(path){
  factors <- colnames(BDLdata)
  Nf = length(factors)
  for (k in 1:Nf){
    # cat(sprintf("%s / %s\n", k, Nf))
    name <- factors[k]
    plot_single_factor(name, path, k)
  }
}

#' Plot of a single factor.
#' 
#' Plots the mean time course and individual samples.
#' @export
plot_single <- function(name_A){
  # data and mean data
  dA <- BDLdata[[name_A]]
  dA.mean <- BDLmean[[name_A]]
  # title
  title <- sprintf("%s (%s)", name_A, BDLfactors$ftype[BDLfactors$id==name_A])
  # basic plot
  plot(BDLsamples$time_fac, dA, xlab="time", ylab=name_A, main=title, col=rgb(0.5,0.5,0.5, 0.4),
       ylim=c(0, max(dA, na.rm=TRUE)*1.1))
  points(BDLsamples$time_fac, dA, col="black")
  points(BDLsamples$time_fac, dA, col=rgb(0,0,1,0.6), pch=16)
  # plot the repeat number
  dA[is.na(dA)] <- -1
  require(calibrate)
  textxy(BDLsamples$time_point, dA, BDLsamples$sid, col="black", cex=0.6)
  
  points(1:nrow(BDLmean), dA.mean, col="red", pch=15)
  lines(1:nrow(BDLmean), dA.mean, col="red")
}


#' Correlation plot between two factors.
#' 
#' Plot the two individual factors and the additional correlation plot
#' for the data.
#' @export
plot_cor_pair <- function(name_A, name_B, single_plots=TRUE){
  require(calibrate)
  if (single_plots){
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
  }
  # correlation plot
  dA <- BDLdata[[name_A]]
  dB <- BDLdata[[name_B]]
  # dA <- dA[!is.na(dA)]
  # d <- dA[!is.na(dA)]
  dA.mean <- BDLmean[[name_A]]
  dB.mean <- BDLmean[[name_B]]
  
  value.max <- 1.1* c(max(dA), max(dB))
  plot(dA, dB, xlim=c(0, value.max[1]), ylim=c(0, value.max[2]), 
       main=sprintf("%s ~ %s", name_A, name_B),
       xlab=name_A, ylab=name_B)
  textxy(dA, dB, BDLsamples$time)
  points(dA.mean, dB.mean, pch=16, col=rgb(0,0,1,0.8) )
  lines(dA.mean, dB.mean, col=rgb(0,0,1,0.8) )
  
  textxy(dA.mean, dB.mean, BDLmean.time, col="blue")
  abline(a=0, b=mean(dB.mean)/mean(dA.mean), col="darkgray")
  
  text(x=0.05*value.max[1], y=0.9*value.max[2], cex=0.5,
       labels=sprintf("S: %1.3f\nS mean: %1.3f\nP: %1.3f\nP mean: %1.3f", 
                      cor(dA, dB, method="spearman", use="pairwise.complete.obs"),
                      cor(dA.mean, dB.mean, method="spearman", use="pairwise.complete.obs"),
                      cor(dA, dB, method="pearson", use="pairwise.complete.obs"),
                      cor(dA.mean, dB.mean, method="pearson", use="pairwise.complete.obs")))
  
  # single plots
  if (single_plots){
    plot_single(name_A)
    plot_single(name_B)
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

#' Heatmap color ramp function.
#' 
#' Definition of heatmap colors via a colorRampPalette.
#' The returned function can be called with the numbers of colors in the heatmap.
#' @export
HeatmapColors <- function(){
  cols <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                     "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")) 
  return(cols)
} 

