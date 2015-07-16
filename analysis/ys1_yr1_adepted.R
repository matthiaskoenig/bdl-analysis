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
}

# Calculate the correlation of the matrix with multiple repeats.
# Works in case of the timecourse and timepoint repeat data.
son.fS_star.matrix <- function(a,b){
  # create vector of feature matrix
  a_vec <- as.vector(a)
  b_vec <- as.vector(b)
  r = son.fS_star(a_vec, b_vec)
}

# Adapted Pearson correlation R* (Shifted to the interval [0, 1])
son.fR_star <- function(a,b){
  r = (cor(x=a, y=b, method='pearson') + 1) / 2
}
# fR_star for matrix data (timecourse and time point repeats)
son.fR_star_matrix <- function(a,b){
  r = (cor(x=as.vector(a), y=as.vector(b), method='pearson') + 1) / 2
}

# --------------------------------------------------------------------------------------------
# Slope based component (adjacent timepoints)
# --------------------------------------------------------------------------------------------

# Slope between two adjacent measurments
son.slope <- function(x1,x2,t1,t2){
  s = ((x2 - x1) / (t2 - t1))  
}

# Sign of slope
son.L <- function(x1,x2,t1,t2)
  L = sign(son.slope(x1,x2,t1,t2));
end

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
    for (q in 1:Nr)
      Aab = Aab + son.fA(a[,p], b[,q], time_pts) 
    }
  }
  r = Aab/Nr/Nr
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

# --------------------------------------------------------------------------------------------
# YS1
# --------------------------------------------------------------------------------------------

# ys1
# the time course vectors for factor a and b and the time_pts are required,
# with all 3 vectors having the same length.
ys1 = function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
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

# 
ys1.timecourse = function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
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


# Calculate ys1 matrix for a given data frame.
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
      # cat(sprintf("[%s, %s]\n", k, i))
      res <- ys1(data[,k], data[,i], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A.mat[k,i] = res$A
      M.mat[k,i] = res$M
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, A=A.mat, M=M.mat))
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
# -----------------------------------------------------------------------

# ----------------------------------------

time
samples$time

name_A = "Ppara"
name_B = "Cyp3a11"
f_cor_pair_plot(name_A, name_B)
ys1(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)


