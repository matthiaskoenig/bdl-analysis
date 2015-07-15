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
##############################################################################################


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
  r = (cor(x=a,y=b, method='spearman') + 1) / 2  
}

# Adapted Pearson correlation R* (Shifted to the interval [0, 1])
son.fR_star <- function(a,b){
  r = (cor(x=a, y=b, method='pearson') + 1) / 2
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


time
samples$time

name_A = "Ppara"
name_B = "Cyp3a11"
f_cor_pair_plot(name_A, name_B)
ys1(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)


