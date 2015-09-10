##############################################################################################
#
#    YS and YR correlation measures
#
#    Correlation score for time course data based on 
#    "A modfied correlation coefficient based similarity measure for clustering
#     time-course gene expression data" [Son2007]
#
#    For recent reviews of correlation scores for gene expression data look
#    Ja
#
#     ys1 and yr1 are a weighted score including contributions of classical correlation score and
#     shape of the time course profile, namely the slope between measurements and the minimal
#     and maximal values."
#
#     Two main use cases are not covered by YS1 and YR1.
#     A] The existance of muliple repeats per time course, i.e. for one sample all time points
#        where measured.
#     B] The existance of multiple repeats per time point, i.e. every time point belongs to a 
#        different sample with time points measured multiple times. 
#     For the analysis of experimental data of class 2 in mouse liver we provide an extended
#     YS1 and YR1 approach applicable to experimental data of class A and B.
#
##############################################################################################

# --------------------------------------------------------------------------------------------
# Classical correlation component (Spearman/Pearson)
# --------------------------------------------------------------------------------------------

#' Adapted Spearman correlation S*.
#' 
#' Shifted spearman correlation on the individual data shifted
#' to the interval [0, 1].
#' @export
son.fS_star <- function(a,b){
  r <- (cor(x=a,y=b, method='spearman') + 1) / 2 
  return(r)
}

#' Adapted Pearson correlation R*.
#' 
#' Shifted pearson correlation on individual data shifted 
#' to the interval [0, 1].
#' @export
son.fR_star <- function(a,b){
  r = (cor(x=a, y=b, method='pearson') + 1) / 2
  return(r)
}

# --------------------------------------------------------------------------------------------
# Slope based component (adjacent timepoints)
# --------------------------------------------------------------------------------------------

#' Slope between two adjacent measurments.
#' 
#' @export
son.slope <- function(x1,x2,t1,t2){
  s = ((x2 - x1) / (t2 - t1))
  return(s)
}

#' Sign of slope.
#' 
#' @export
son.L <- function(x1,x2,t1,t2){
  L = sign(son.slope(x1,x2,t1,t2))
  return(L)
}

#' Relative count of number of identical slopes between two time courses.
#' 
#' @export
son.fA <- function(a,b,time_pts){
  N = length(a)
  count = 0
  for (i in 1:(N-1)){
    if (son.L(a[i],a[i+1],time_pts[i],time_pts[i+1]) == son.L(b[i],b[i+1],time_pts[i],time_pts[i+1])){
      count = count + 1
    }
  }
  r = count / (N-1)
  return(r)
}

#' Adapted slope measurement.
#'
#' Since the concordance index A counts only the number of time intervals
#' in which there is an agreement in the sign of the change between two
#' profiles, it may loose the information of the size of the change.
#' Use the size of the change and compute the correlation coefficient between the
#' differences of the profiles in the time intervals.
#' This works on the distances, not the slopes! Which is different if the time points
#' are not equidistant.
#' @export
son.fA_star <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)]
  dB = b[2:Nt]-b[1:(Nt-1)]
  r = (cor(x=dA, y=dB, method='pearson') + 1) / 2
  return(r)
}

#' Adapted slope measurement (spearman)
#' 
#' Similar to fA_star but uses spearman for slope correlation calculation.
#' @export
son.fA_star2 <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)] 
  dB = b[2:Nt]-b[1:(Nt-1)]
  r = (cor(x=dA, y=dB, method='spearman') + 1) / 2
  return(r)
}

#' Calculating the distance on the slopes.
#' 
#' The distance based measurement should be preferred in case of unequal time points.
#' @export
son.fA_star3 <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)] 
  dB = b[2:Nt]-b[1:(Nt-1)]
  dt = time_pts[2:Nt]-time_pts[1:(Nt-1)]
  r = (cor(x=dA/dt, y=dB/dt, method='pearson') + 1) / 2
  return(r)
}

# --------------------------------------------------------------------------------------------
# Min/Max component (based on global time course information)
# --------------------------------------------------------------------------------------------

#' Compare timepoint of minimal and maximal value.
#' 
#' @export
son.fM = function(a,b){
  # find index of max and min
  idx_max_a = which.max(a)
  idx_max_b = which.max(b)
  idx_min_a = which.min(a)
  idx_min_b = which.min(b)
  
  # compare max and min indices
  if ( (idx_max_a == idx_max_b) && (idx_min_a == idx_min_b) ){
    r = 1.0
  } else if ( (idx_max_a == idx_max_b) || (idx_min_a == idx_min_b) ){
    r = 0.5
  } else if ( (idx_max_a != idx_max_b) && (idx_min_a != idx_min_b) ){
    r = 0.0
  }
  return(r)
}

#' Alternative index to M.
#'
#' M index utilizing the distances between two profiles time 
#' points where the max/min is attained.
#' Uses the difference in indices !
#' @export
son.fM_star = function(a,b,time_pts){
  Nt = length(time_pts)
  idx_max_a = which.max(a)[1]
  idx_max_b = which.max(b)[1]
  idx_min_a = which.min(a)[1]
  idx_min_b = which.min(b)[1]
  r = 1 - ( (abs(idx_min_a-idx_min_b) + abs(idx_max_a-idx_max_b))/(2*(Nt-1))  )
  return(r)
}

#' Alternative index to M.
#' 
#' Alternative index to M which utilizes the distances between two profiles time.
#' points where the max/min is attained.
#' Uses the difference in times
#' @export
son.fM_star2 = function(a,b,time_pts){
  Nt = length(time_pts)
  idx_max_a = which.max(a)[1]
  idx_max_b = which.max(b)[1]
  idx_min_a = which.min(a)[1]
  idx_min_b = which.min(b)[1]
  r = 1 - ( (abs(time_pts[idx_min_a]-time_pts[idx_min_b]) + abs(time_pts[idx_max_a]-time_pts[idx_max_b]))/(2*(time_pts[Nt]-time_pts[1])))
  return(r)
}

#' Curvature based M measure. 
#'
#' Contains implicitly the information about minimum and maximum 
#' of a timecourse.
#' @export 
son.fM_star3 = function(a,b,time_pts){
  Nt = length(time_pts)
  # first derivative (slope)
  dA = a[2:Nt]-a[1:(Nt-1)] 
  dB = b[2:Nt]-b[1:(Nt-1)]
  dt = time_pts[2:Nt]-time_pts[1:(Nt-1)]
  # second derivative (curvature)
  d2A = dA[2:(Nt-1)]-dA[1:(Nt-2)] 
  d2B = dB[2:(Nt-1)]-dB[1:(Nt-2)]
  # mid values
  dtm = (time_pts[2:Nt]+time_pts[1:(Nt-1)])/2
  d2t = dtm[2:(Nt-1)]-dt[1:(Nt-2)]
  r = (cor(x=d2A/d2t, y=d2B/d2t, method='pearson') + 1) / 2
  return(r)
}

# --------------------------------------------------------------------------------------------
# Calculation of all component matrices
# --------------------------------------------------------------------------------------------

#' YS and YR correlation measure.
#' 
#' The time course vectors for factor a and b and the time_pts are required,
#' with all 3 vectors having the same length.
#' This function only calculates the components of the various correlation 
#' scores which still have to be combined via respective weights for the
#' components, i.e. for instance
#' w1=0.50, w2=0.25, w3=0.25
#' @export
ysr <- function(a,b,time_pts, use="all.obs"){
  # check if any time point problems
  if (any(table(time_pts)>1)){
    stop("time_pts not unique")
  }  
  # create a, b data frame
  df <- data.frame(a, b, time_pts)
  
  # handle NAs
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }  
  if (identical(use, "pairwise.complete.obs")){
    df <- df[complete.cases(df),]
  }  
  # correlation
  R_star <- son.fR_star(df$a,df$b)
  S_star <- son.fS_star(df$a,df$b)
  # slope
  A <- son.fA(df$a,df$b,df$time_pts)
  A_star <- son.fA_star(df$a,df$b,df$time_pts)
  A_star2 <- son.fA_star2(df$a,df$b,df$time_pts)
  # min/max
  M <- son.fM(df$a,df$b)
  M_star <- son.fM_star(df$a,df$b, df$time_pts)
  M_star2 <- son.fM_star2(df$a,df$b, df$time_pts)
  M_star3 <- son.fM_star3(df$a,df$b, df$time_pts)
  
  return( list(R_star=R_star,
               S_star=S_star, 
               A=A, 
               A_star=A_star,
               A_star2=A_star2,
               M=M,
               M_star=M_star,
               M_star2=M_star2,
               M_star3=M_star3) )
}


#' Calculate ys and yr component matrices given data frame.
#' 
#' Every column of the data.frame is a single measurement (or mean) of a single factor for
#' the provided time points. 
#' Use this function to calculate the correlation matrix on the mean data.
#' @export
ysr.matrices <- function(data, time_pts, use="pairwise.complete.obs"){
  N <- ncol(data)
  value.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(value.mat) <- names(data)
  rownames(value.mat) <- names(data)
  
  # all components used for the various scors
  R_star <- value.mat
  S_star <- value.mat
  A <- value.mat
  A_star <- value.mat
  A_star2 <- value.mat
  M <- value.mat
  M_star <- value.mat
  M_star2 <- value.mat
  M_star3 <- value.mat
  
  for (k in 1:N){
    for (i in 1:N){
      # calculate the simple score pairwise between all factors
      res <- ysr(data[,k], data[,i], time_pts, use=use)
      
      # correlation
      R_star[k,i] = res$R_star
      S_star[k,i] = res$S_star
      # slope
      A[k,i] = res$A
      A_star[k,i] = res$A_star
      A_star2[k,i] = res$A_star2
      # min/max
      M[k,i] = res$M
      M_star[k,i] = res$M_star
      M_star2[k,i] = res$M_star2
      M_star3[k,i] = res$M_star3
    }
  }
  return(list(R_star=R_star,
              S_star=S_star, 
              A=A,
              A_star=A_star,
              A_star2=A_star2,
              M=M,
              M_star=M_star,
              M_star2=M_star2,
              M_star3=M_star3))
}
