# ---- ys1 and yr1 definition -------------------------------------------------------------------
# "A modfied correlation coefficient based similarity measure for clustering
#  time-course gene expression data" {Son2007}
# ys1 and yr1 are a weighted score including contributions of classical correlation score and
# shape of the time course profile, namely the slope between measurements and the minimal
# and maximal values.

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

# Adapted Spearman correlation S* (Shifted to the interval [0, 1])
son.fS_star <- function(a,b){
  r = (cor(x=x1,y=x2, method='spearman') + 1) / 2  
}

# Adapted Pearson correlation R* (Shifted to the interval [0, 1])
son.fR_star <- function(a,b){
  r = (cor(x=x1, y=x2, method='pearson') + 1) / 2
}

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

# ys1
# the time course vectors for factor a and b and the time_pts are required,
# with all 3 vectors having the same length.
ys1 = function(a,b,time_pts, w1=0.50, w2=0.25, w3=0.25, use="all.obs"){
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)){
    stop("invalid 'use' argument")
  }
  df <- data.frame(a, b, time_pts)
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

# -----------------------------------------------------------------------
# implementation tests from [Son2008]
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

# -----------------------------------------------------------------------