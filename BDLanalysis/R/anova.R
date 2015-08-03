


#' Perform anova for given factor.
#' 
#' Similar implementation to
#' http://www.r-tutor.com/elementary-statistics/analysis-variance/completely-randomized-design
#' @export
single_factor_anova <- function(name){
  # data matrix
  mat1 <- t(BDLmatrices[[name]])
  colnames(mat1) <- levels(BDLsamples$time_fac)
  
  # Concatenate the data rows of df1 into a single vector r .
  r = c(t(as.matrix(mat1))) # response data 
  
  # Assign new variables for the treatment levels and number of observations.
  f = levels(BDLsamples$time_fac)   # treatment levels 
  k = 8                          # number of treatment levels 
  n = 5                          # observations per treatment 
  
  # Create a vector of treatment factors that corresponds to each element of r in step 3 with the gl function.
  tm = gl(k, 1, n*k, factor(f))   # matching treatments 
  
  # Apply the function aov to a formula that describes the response r by the treatment factor tm.
  # Fit an analysis of variance model
  av = aov(r ~ tm) 
  
  # Print out the ANOVA table with the summary function. 
  # summary(av)
  return(av)
}

#' Performs ANOVA for all factors.
#' 
#' ANOVA for all single factors are calculated. 
#' The reported p-values are tghe unadjusted p-values (no multiple testing
#' correction)
#' @export
all_factor_anova <- function(){
  factors <- colnames(BDLdata)
  df.anova <- data.frame(factors)
  df.anova$p.value <- NA
  
  for (k in 1:nrow(df.anova)){
    name <- factors[[k]]
    av <- single_factor_anova(name)
    p.value <- summary(av)[[1]][["Pr(>F)"]][[1]]
    df.anova[k, "p.value"] <- p.value
  }
  return(df.anova)
}

#' Create significant code for given p.value.
#' 
#' Function returns a significance string for the 
#' given p-value similar to the internal R functions.
#' @export
significant_code <- function(p.value){
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  sig = " "
  if (p.value <= 0.001){
    sig = "***"
  } else if (p.value <= 0.01){
    sig = "**"
  } else if (p.value <=0.05){
    sig = "*"
  } else if (p.value <=0.1){
    sig = "."
  } else if (p.value <=1){
    sig = " "
  }
  return (sig)
}