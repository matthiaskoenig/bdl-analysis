##############################################################################################
#
#    YS1 and YR1 adapted
#
#    Matthias Koenig
#    2015-07-13
#
#    YS1 and YR1 implementation based on
#   "A modfied correlation coefficient based similarity measure for clustering
#     time-course gene expression data" [Son2007]
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

# Adapted Spearman correlation S* (Shifted to the interval [0, 1])
son.fS_star <- function(a,b){
  r = (cor(x=a,y=b, method='spearman') + 1) / 2 
  return(r)
}

# Calculate the correlation of the matrix with multiple repeats.
# Works in case of the timecourse and timepoint repeat data.
son.fS_star.matrix <- function(a,b){
  # create vector of feature matrix
  a_vec <- as.vector(a)
  b_vec <- as.vector(b)
  r = son.fS_star(a_vec, b_vec)
  return(r)
}

# Adapted Pearson correlation R* (Shifted to the interval [0, 1])
son.fR_star <- function(a,b){
  r = (cor(x=a, y=b, method='pearson') + 1) / 2
  return(r)
}
# fR_star for matrix data (timecourse and time point repeats)
son.fR_star_matrix <- function(a,b){
  r = (cor(x=as.vector(a), y=as.vector(b), method='pearson') + 1) / 2
  return(r)
}

# --------------------------------------------------------------------------------------------
# Slope based component (adjacent timepoints)
# --------------------------------------------------------------------------------------------

# Slope between two adjacent measurments
son.slope <- function(x1,x2,t1,t2){
  s = ((x2 - x1) / (t2 - t1))
  return(s)
}

# Sign of slope
son.L <- function(x1,x2,t1,t2){
  L = sign(son.slope(x1,x2,t1,t2))
  return(L)
}

# Relative count of number of identical slopes between two time courses.
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

# Since the concordance index A counts only the number of time intervals
# in which there is an agreement in the sign of the change between two
# profiles, it may loose the information of the size of the change.
# Use the size of the change and compute the correlation coefficient between the
# differences of the profiles in the time intervals.
# This works on the distances, not the slopes! Which is different if the time points
# are not equidistant.
son.fA_star <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)]
  dB = b[2:Nt]-b[1:(Nt-1)]
  r = (cor(x=dA, y=dB, method='pearson') + 1) / 2
  return(r)
}

# Calculating the distance on the slopes, which should be 
# preferred in case of unequal time points
son.fA_star2 <- function(a,b,time_pts){
  Nt = length(time_pts)
  dA = a[2:Nt]-a[1:(Nt-1)] 
  dB = b[2:Nt]-b[1:(Nt-1)]
  dt = time_pts[2:Nt]-time_pts[1:(Nt-1)]
  r = (cor(x=dA/dt, y=dB/dt, method='pearson') + 1) / 2
  return(r)
}


# fA calculation for time courses with muliple repeats,
# i.e. multiple samples measured for the same time course.
# a and b are matrices of the form [Nt, Nr], with
#   Nt : number of time points
#   Nr : number of repeats == number of samples per condition in case of time course
# This is only applicable if the individual measurements within a column are from the
# same sample !, i.e. real timecourse repeats. If every timepoint within a repeat 
# belongs to a different sample the timepoint version must be used.
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

# Calculate all the permutations of the datapoint slopes,
# and count the identical ones.
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

# Compare timepoint of minimal and maximal value.
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

# Alternative index to M which utilizes the distances between two profiles time 
# points where the max/min is attained.
# Uses the difference in indices !
son.fM_star = function(a,b,time_pts){
  Nt = length(time_pts)
  idx_max_a = which.max(a)[1]
  idx_max_b = which.max(b)[1]
  idx_min_a = which.min(a)[1]
  idx_min_b = which.min(b)[1]
  r = 1 - ( (abs(idx_min_a-idx_min_b) + abs(idx_max_a-idx_max_b))/(2*(Nt-1))  )
  return(r)
}

# Using the curvature which contains the information about min/max but is more
# robust 
son.fM_star2 = function(a,b,time_pts){
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
  # print(d2A/d2t)
  # print("")
  # print(d2B/d2t)
  r = (cor(x=d2A/d2t, y=d2B/d2t, method='pearson') + 1) / 2
  return(r)
}


# fM calculation for 
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

# ys1
# the time course vectors for factor a and b and the time_pts are required,
# with all 3 vectors having the same length.
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

# ys1 for data matrices with repeated time course measurements in columns.
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

# ys1 calculation for timepoint data. Every measurement belongs to one sample.
# The columns are the number of repeats.
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

# --------------------------------------------------------------------------------------------
# YS2
# --------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------
# YR1
# --------------------------------------------------------------------------------------------

# yr1
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

# Calculate yr1 matrix for given data frame.
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

# -----------------------------------------------------------------------
# Implementation Tests from [Son2008]
# -----------------------------------------------------------------------
time_pts <- seq(0,4)
x1 <- c(2,3,6,4,7)
x2 <- c(1,2,3,5,3)
y1 <- c(4,3,6,2,7)
y2 <- c(5,2,3,1,3)

# plot data
par(mfrow=c(1,2))
plot(time_pts, x1, col="blue", pch=16, type="o", ylim=c(0, max(max(x1), max(x2))),
     main = "x1, x2", xlab="values", ylab="time")
points(time_pts, x2, col="red", pch=15, type="o")

plot(time_pts, y1, col="blue", pch=16, type="o", ylim=c(0, max(max(y1), max(y2))),
     main = "y1, y2", xlab="values", ylab="time")
points(time_pts, y2, col="red", pch=15, type="o")
par(mfrow=c(1,1))

# ys1
test_res <- function(y, y_expected, digits=3){
  stopifnot(round(y, digits=digits) == y_expected)
}

# correlation values
test_res(cor(x=x1,y=x2, method='spearman'), 0.667)
test_res(cor(x=x1,y=x2, method='pearson'),  0.439)
# ys1
test_res(ys1(x1, x2, time_pts, w1=0.25, w2=0.50, w3=0.25)$value, 0.583)
test_res(ys1(y1, y2, time_pts, w1=0.25, w2=0.50, w3=0.25)$value, 0.833)
# yr1
test_res(yr1(x1, x2, time_pts, w1=0.25, w2=0.50, w3=0.25)$value, 0.555)
test_res(yr1(y1, y2, time_pts, w1=0.25, w2=0.50, w3=0.25)$value, 0.805)

rm(test_res, x1, x2, y1, y2, time_pts)

# TODO: implement the tests for ys2
# -----------------------------------------------------------------------



