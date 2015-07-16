##############################################################################################
#
#    YS1 and YR1 adapted
#
#    Matthias Koenig
#    2015-07-13
#
#    Helper functions to apply ys1 and yr1 to the BDL dataset.
#
##############################################################################################

source("ys1_yr1.R")  # definition of ys1 and yr1

# Calculate ys1 matrix for a given data frame.
# Every column of the data.frame is a single measurement (or mean) of a single factor for
# the provided time points. 
# Use this function to calculate the correlation matrix on the mean data.
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

ys2.df <- function(data, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
  N <- ncol(data)
  value.mat <- matrix(NA, nrow=N, ncol=N)
  colnames(value.mat) <- names(data)
  rownames(value.mat) <- names(data)
  S_star.mat <- value.mat
  A_star.mat <- value.mat
  M_star.mat <- value.mat
  
  for (k in 1:N){
    for (i in 1:N){
      # calculate the simple score pairwise between all factors
      res <- ys2(data[,k], data[,i], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A_star.mat[k,i] = res$A_star
      M_star.mat[k,i] = res$M_star
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, A_star=A_star.mat, M=M_star.mat))
}



# Calculate ys1 for factors provided as list of matrices with size [Nt, Nr].
ys1.timecourse.df <- function(data_list, time_pts, w1=0.50, w2=0.25, w3=0.25, use="pairwise.complete.obs"){
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
      res <- ys1.timecourse(data_list[[k]], data_list[[i]], time_pts, w1=w1, w2=w2, w3=w3, use=use)
      value.mat[k,i] = res$value
      S_star.mat[k,i] = res$S_star
      A.mat[k,i] = res$A
      M.mat[k,i] = res$M
    }
  }
  return(list(value=value.mat, S_star=S_star.mat, A=A.mat, M=M.mat))
}


# Calculate ys1 for factors provided as list of matrices with size [Nt, Nr].
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

ys1.timepoints(data_list[[1]], data_list[[2]], time_pts=dmean.time)
ys1.timepoints(data_list[[136]], data_list[[42]], time_pts=dmean.time)

ys1.timecourse(data_list[[1]], data_list[[2]], time_pts=dmean.time)
ys1.timecourse.df(data_list[1:2], time_pts=dmean.time)

data_list[1:2]
name_A = "Ppara"   # [[1]]
name_B = "Cyp3a11" # [[2]]
f_cor_pair_plot(name_A, name_B)
ys1(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)
ys2(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)

# 136 Actb
# 42 Actb.x
# 89 Actb.y

ys1.timecourse.df(data_list[c(42,89,136)], time_pts=dmean.time)

ys1.timecourse.df(data_list[c(1,154)], time_pts=dmean.time)


#--------------------------------------------------------------------------------------------------------
# Example calculations on BDL dataset.

# calculate ys1 on mean timecourse data
ys1.res <- ys1.df(dmean, dmean.time)

source("ys1_yr1.R")  # definition of ys1 and yr1

ys2.res <- ys2.df(dmean, dmean.time)
f_corrplot("cor.ys2", data=ys2.res$value, order="original")
f_corrplot("cor.ys2-scaled", data=2*(ys2.res$value-0.5), order="original")
f_corrplot("cor.ys2", data=ys2.res$value, order="hclust")
f_corrplot("cor.ys2-scaled", data=2*(ys2.res$value-0.5), order="hclust")
options$height=1600
options$width=1600
summary(ys2.res$value)

