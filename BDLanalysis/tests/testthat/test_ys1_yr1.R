##############################################################################################
#    Implementation tests for YS and YR
##############################################################################################

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
