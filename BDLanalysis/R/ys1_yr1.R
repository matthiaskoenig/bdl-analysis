##############################################################################################
#
#    YS1 and YR1 correlation measures
#
#    Matthias Koenig
#    2015-08-03
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

# TODO: implement how to handle small fluctuations in slope.

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


#' fS_star calculation for matrix with multiple repeats.
#' 
#' Works in case of the timecourse and timepoint repeat data.
#' @export
son.fS_star.matrix <- function(a,b){
  # create vector of feature matrix
  a_vec <- as.vector(a)
  b_vec <- as.vector(b)
  r = son.fS_star(a_vec, b_vec)
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


#' fR_star calculation for matrix with multiple repeats.
#' 
#' Works for timecourse and timepoint repeat data.
#'@export
son.fR_star.matrix <- function(a,b){
  r = (cor(x=as.vector(a), y=as.vector(b), method='pearson') + 1) / 2
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


#' Calculating the distance on the slopes.
#' 
#' The distance based measurement should be preferred in case of unequal time points.
#' @export
son.fA_star2 <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)] 
  dB = b[2:Nt]-b[1:(Nt-1)]
  dt = time_pts[2:Nt]-time_pts[1:(Nt-1)]
  r = (cor(x=dA/dt, y=dB/dt, method='pearson') + 1) / 2
  return(r)
}


#' fA calculation for time courses with multiple repeats.
#' 
#' fA with multiple repeats i.e. multiple samples measured for the same time course.
#' a and b are matrices of the form [Nt, Nr], with
#'   Nt : number of time points
#'   Nr : number of repeats == number of samples per condition in case of time course
#' This is only applicable if the individual measurements within a column are from the
#' same sample !, i.e. real timecourse repeats. If every timepoint within a repeat 
#' belongs to a different sample the timepoint version must be used.
#' @export
son.fA.timecourse <- function(a,b,time_pts){
  Nr = ncol(a)
  Aab = 0
  # calulate the mean A over all repeats of time course
  for (p in 1:Nr){
    for (q in 1:Nr){
      Aab = Aab + son.fA(a[,p], b[,q], time_pts) 
    }
  }
  r = Aab/Nr/Nr
  return(r)
}


#' fA caluclation for time points with multiple repeats.
#' 
#' Calculate all the permutations of the datapoint slopes,
#' and count the identical ones.
#' @export
son.fA.timepoints <- function(a,b,time_pts){
  Nt = length(time_pts)
  Nr = ncol(a)
  count = 0
  for (i in 1:(Nt-1)){
    for (p in 1:Nr){
      for (q in 1:Nr){
        if (son.L(a[i,p],a[i+1,p],time_pts[i],time_pts[i+1]) == son.L(b[i,q],b[i+1,q],time_pts[i],time_pts[i+1])){
          count = count + 1
        }
      }
    }
  }
  r = count/(Nt-1)/Nr/Nr
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


#' fM calculation for timecourse.
#' 
#' @export
son.fM.timecourse <- function(a,b){
  Nr = ncol(a)
  Mab = 0
  # calulate the mean A over all repeats of time course
  for (p in 1:Nr){
    for (q in 1:Nr){
      Mab = Mab + son.fM(a[,p], b[,q]) 
    }
  }
  r = Mab/Nr/Nr
  return(r)
}

# --------------------------------------------------------------------------------------------
# YS1
# --------------------------------------------------------------------------------------------

#' ys1 correlation measure.
#' 
#' the time course vectors for factor a and b and the time_pts are required,
#' with all 3 vectors having the same length.
#' @export
ys1 <- function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
  
  # check if any
  if (any(table(time_pts)>1)){
    stop("time_pts not unique")
  }  
  # handle NAs
  if (identical(use, "pairwise.complete.obs")){
    df <- df[complete.cases(df),]
  }  
  S_star <- son.fS_star(df$a,df$b)
  A <- son.fA(df$a,df$b,df$time_pts)
  M <- son.fM(df$a,df$b)
  value = w1*S_star + w2*A + w3*M
  return( list(value=value, S_star=S_star, A=A, M=M) )
}


#' ys1 for time course data.
#' 
#' Uses data matrices with repeated time course measurements in columns.
#' @export
ys1.timecourse = function(a,b, time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
  # check if any
  if (any(table(time_pts)>1)){
    stop("time_pts not unique")
  }  
  # handle NAs
  # TODO: handle the complete cases in matrix
  # if (identical(use, "pairwise.complete.obs")){
  #   df <- df[complete.cases(df),]
  # }  
  
  S_star <- son.fS_star.matrix(a, b)
  A <- son.fA.timecourse(a, b, time_pts)
  M <- son.fM.timecourse(a, b)
  value = w1*S_star + w2*A + w3*M
  return( list(value=value, S_star=S_star, A=A, M=M) )
}


#' ys1 calculation for timepoint data. 
#' 
#' Every measurement belongs to one sample. The columns are the number of repeats.
#' @export
ys1.timepoints = function(a,b,time_pts, w1=0.50, w2=0.5, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
  # check if any
  if (any(table(time_pts)>1)){
    stop("time_pts not unique")
  }
  # handle NAs
  # TODO: handle the complete cases in matrix
  # if (identical(use, "pairwise.complete.obs")){
  #   df <- df[complete.cases(df),]
  # }  
  
  S_star <- son.fS_star.matrix(a, b)
  A <- son.fA.timepoints(a, b, time_pts)
  value = w1*S_star + w2*A
  return( list(value=value, S_star=S_star, A=A) )
}


#' Calculate ys1 matrix for a given data frame.
#' 
#' Every column of the data.frame is a single measurement (or mean) of a single factor for
#' the provided time points. 
#' Use this function to calculate the correlation matrix on the mean data.
#' TODO: refactor so that only one function for the calculation of correlation matrix exists.
#' @export
ys1.df <- function(data, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  N <- ncol(data)
  value.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(value.mat) <- names(data)
  rownames(value.mat) <- names(data)
  S_star.mat <- value.mat
  A.mat <- value.mat
  M.mat <- value.mat
  
  for (k in 1:N){
    for (i in 1:N){
      # calculate the simple score pairwise between all factors
      res <- ys1(data[,k], data[,i], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A.mat[k,i] = res$A
      M.mat[k,i] = res$M
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, A=A.mat, M=M.mat))
}

#' Calculate ys1 for factors provided as list of matrices with size [Nt, Nr].
#' 
#' @export
ys1.timepoints.df <- function(data_list, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  Nf <- length(data_list)
  value.mat <- matrix(NA, nrow=Nf, ncol=Nf)
  colnames(value.mat) <- names(data_list)
  rownames(value.mat) <- names(data_list)
  S_star.mat <- value.mat
  A.mat <- value.mat
  M.mat <- value.mat
  
  for (k in 1:Nf){
    for (i in 1:Nf){
      cat(sprintf("[%s, %s]\n", k, i))
      res <- ys1.timepoints(data_list[[k]], data_list[[i]], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A.mat[k,i] = res$A
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, A=A.mat))
}

# --------------------------------------------------------------------------------------------
# YS2
# --------------------------------------------------------------------------------------------

#' ys2 calculation.
#' 
#' @export
ys2 <- function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
  
  # check if any
  if (any(table(time_pts)>1)){
    stop("time_pts not unique")
  }  
  # handle NAs
  if (identical(use, "pairwise.complete.obs")){
    df <- df[complete.cases(df),]
  }  
  S_star <- son.fS_star(df$a, df$b)
  A_star <- son.fA_star(df$a, df$b, df$time_pts)
  A_star2 <- son.fA_star2(df$a, df$b, df$time_pts)
  M_star <- son.fM_star(df$a, df$b, df$time_pts)
  M_star2 <- son.fM_star2(df$a, df$b, df$time_pts)
  value = w1*S_star + w2*A_star + w3*M_star
  return( list(value=value, S_star=S_star, 
               A_star=A_star, A_star2=A_star2, 
               M_star=M_star, M_star2=M_star2) )
}


#' ys2 correlation matrix.
#' 
#' TODO: refactor in one function
#' 
#' @export
ys2.df <- function(data, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  N <- ncol(data)
  value.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(value.mat) <- names(data)
  rownames(value.mat) <- names(data)
  S_star.mat <- value.mat
  A_star.mat <- value.mat
  A_star2.mat <- value.mat
  M_star.mat <- value.mat
  M_star2.mat <- value.mat
  for (k in 1:N){
    for (i in 1:N){
      # calculate the simple score pairwise between all factors
      res <- ys2(data[,k], data[,i], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A_star.mat[k,i] = res$A_star
      A_star2.mat[k,i] = res$A_star2
      M_star.mat[k,i] = res$M_star
      M_star2.mat[k,i] = res$M_star2
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, 
              A_star=A_star.mat, A_star2=A_star2.mat, 
              M_star=M_star.mat, M_star2=M_star2.mat))
}


# --------------------------------------------------------------------------------------------
# YR1
# --------------------------------------------------------------------------------------------

#' yr1 calculation.
#' 
#' @export
yr1 = function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
  # handle NAs
  if (identical(use, "pairwise.complete.obs")){
    df <- df[complete.cases(df),]
  }  
  R_star <- son.fR_star(df$a,df$b)
  A <- son.fA(df$a,df$b,df$time_pts)
  M <- son.fM(df$a,df$b)
  value = w1*R_star + w2*A + w3*M
  return( list(value=value, R_star=R_star, A=A, M=M) )
}

#' yr1 correlation matrix for given data frame.
#' 
#' This calculates the yr1 correlation matrix for the given set of factors.
#' @export
yr1.df <- function(data, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  N <- ncol(data)
  cor.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(cor.mat) <- names(data)
  rownames(cor.mat) <- names(data)
  for (k in 1:N){
    for (i in 1:N){
      # cat(sprintf("[%s, %s]\n", k, i))
      cor.mat[k,i] = yr1(data[,k], data[,i], time_pts, w1=w1, w2=w2, w3=w3, use=use)$value
    }
  }
  return(cor.mat)
}
