# -----------------------------------------------------------------------
# Implementation Tests from [Son2008]
# -----------------------------------------------------------------------
time_pts <- seq(0,4)
x1 <- c(2,3,6,4,7)
x2 <- c(1,2,3,5,3)
y1 <- c(4,3,6,2,7)
y2 <- c(5,2,3,1,3)

# # plot data
# par(mfrow=c(1,2))
# plot(time_pts, x1, col="blue", pch=16, type="o", ylim=c(0, max(max(x1), max(x2))),
#      main = "x1, x2", xlab="values", ylab="time")
# points(time_pts, x2, col="red", pch=15, type="o")
# 
# plot(time_pts, y1, col="blue", pch=16, type="o", ylim=c(0, max(max(y1), max(y2))),
#      main = "y1, y2", xlab="values", ylab="time")
# points(time_pts, y2, col="red", pch=15, type="o")
# par(mfrow=c(1,1))

# correlation values
test_that("Correlation between x1 and x2", {
  expect_equal(cor(x=x1,y=x2, method='spearman'), 0.667, tolerance = .001)
  expect_equal(cor(x=x1,y=x2, method='pearson'), 0.439, tolerance = .001)
})

test_that("ys1 with test data", {
  w1=0.25 
  w2=0.50
  w3=0.25
  x.yrs <- ysr(a=x1, b=x2, time_pts=time_pts)
  expect_equal(w1*x.yrs$S_star + w2*x.yrs$A + w3*x.yrs$M, 0.583, tolerance =0.001)
  y.yrs <- ysr(a=y1, b=y2, time_pts=time_pts)
  expect_equal(w1*y.yrs$S_star + w2*y.yrs$A + w3*y.yrs$M, 0.833, tolerance =0.001)
})

test_that("yr1 with test data", {
  w1=0.25 
  w2=0.50
  w3=0.25
  x.yrs <- ysr(a=x1, b=x2, time_pts=time_pts)
  expect_equal(w1*x.yrs$R_star + w2*x.yrs$A + w3*x.yrs$M, 0.555, tolerance=0.001)
  y.yrs <- ysr(a=y1, b=y2, time_pts=time_pts)
  expect_equal(w1*y.yrs$R_star + w2*y.yrs$A + w3*y.yrs$M, 0.805, tolerance =0.001)
})


# TODO: implement the tests for ys2


# 
# ys1.timepoints(data_list[[1]], data_list[[2]], time_pts=dmean.time)
# ys1.timepoints(data_list[[136]], data_list[[42]], time_pts=dmean.time)
# 
# ys1.timecourse(data_list[[1]], data_list[[2]], time_pts=dmean.time)
# ys1.timecourse.df(data_list[1:2], time_pts=dmean.time)
# 
# data_list[1:2]
# name_A = "Ppara"   # [[1]]
# name_B = "Cyp3a11" # [[2]]
# f_cor_pair_plot(name_A, name_B)
# ys1(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)
# ys2(a=dmean[[name_A]], b=dmean[[name_B]], time_pts=dmean.time, w1=0.25, w2=0.5, w3=0.25)
# 
# # 136 Actb
# # 42 Actb.x
# # 89 Actb.y
# 
# ys1.timecourse.df(data_list[c(42,89,136)], time_pts=dmean.time)
# 
# ys1.timecourse.df(data_list[c(1,154)], time_pts=dmean.time)
# 
# #--------------------------------------------------------------------------------------------------------
# # Example calculations on BDL dataset.
# 
# # calculate ys1 on mean timecourse data
# ys1.res <- ys1.df(dmean, dmean.time)
# 
# source("ys1_yr1.R")  # definition of ys1 and yr1
# 
# ys2.res <- ys2.df(dmean, dmean.time)
# f_corrplot("cor.ys2", data=ys2.res$value, order="original")
# f_corrplot("cor.ys2-scaled", data=2*(ys2.res$value-0.5), order="original")
# f_corrplot("cor.ys2", data=ys2.res$value, order="hclust")
# f_corrplot("cor.ys2-scaled", data=2*(ys2.res$value-0.5), order="hclust")
# options$height=1600
# options$width=1600
# summary(ys2.res$value)
# 
