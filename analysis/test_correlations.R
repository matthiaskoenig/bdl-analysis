# Analysis of the correlation structure
n1 <- "Ppara"
n2 <- "Cyp24a1"

ys2.df(data=BDLmean[, c(n1, n2)], time_pts=BDLmean.time)
dev.off()
plot_cor_pair(n1, n2)

Nt <- nrow(BDLmean)
a <- BDLmean[[n1]]
b <- BDLmean[[n2]]

dA = a[2:Nt]-a[1:(Nt-1)]
dB = b[2:Nt]-b[1:(Nt-1)]
dA
dB
plot(dA, dB, col="blue", pch=21)

# A*
r = (cor(x=dA, y=dB, method='pearson') + 1) / 2
s = (cor(x=dA, y=dB, method='spearman') + 1) / 2
s
r

cor(x=dA, y=dB, method='pearson')

# A**
time_pts <- BDLmean.time
dt = time_pts[2:Nt]-time_pts[1:(Nt-1)]
plot(dA/dt, dB/dt)
r = (cor(x=dA/dt, y=dB/dt, method='pearson') + 1) / 2
s = (cor(x=dA/dt, y=dB/dt, method='spearman') + 1) / 2
r
s

# Correlation between distances to min
son.fM_star(a,b,time_pts)



a.norm <- (a-min(a))/(max(a)-min(a))
b.norm <- (b-min(b))/(max(b)-min(b))
plot(a.norm, b.norm)
textxy(a.norm, b.norm, labs = as.character(time_pts), cex=0.8)

# cluster 
cor.cluster[which(groups==1), which(groups==1)]