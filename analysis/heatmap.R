# Example heatmaps with clustering
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/heatmap
source("http://www.bioconductor.org/biocLite.R")
biocLite("ALL")


col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  

heatmap(cor.cluster, col=topo.colors(100))
# already scaled, change colors and size
heatmap(cor.cluster, col=col2(10), scale="none", cexRow=0.5, cexCol=0.5, symm=TRUE)

library("gplots")
heatmap.2(cor.cluster, col=col2(10), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5)

heatmap.2(cor.cluster, col=col2(10), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, density.info="none")
# Display the clusters in addition
heatmap.2(cor.cluster, col=col2(40), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5,
          density.info="none", dendrogram="column", keysize=0.8)

# Perform the clustering
Ngroups = 6
hc <- hclust(dist(cor.cluster)) 
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)

groups


# define colors for 6 cluster (colorbrewer)
library("RColorBrewer")
# display.brewer.all()
colorset <- brewer.pal(Ngroups, "Set1")

color.map <- function(cluster_id) {return(colorset[cluster_id])}
clusterColors <- unlist(lapply(groups, color.map))

png(filename=sprintf("../results/%s", "test.png"), width=options$width, height=options$height, res=options$res)
heatmap.2(cor.cluster, col=col2(100), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5,
          density.info="none", dendrogram="column", Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), keysize=0.8,
          ColSideColors=clusterColors, revC=TRUE)
dev.off()

png(filename=sprintf("../results/%s", "test.png"), width=options$width, height=options$height, res=options$res)
heatmap.2(cor.cluster, col=col2(100), scale="none",
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5,
          density.info="none", dendrogram="column", Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), keysize=0.8,
          ColSideColors=clusterColors, revC=TRUE, colsep=[])
dev.off()
levels(as.factor(groups))


# TODO: add the cluster colors

hc <- hclust(dist(cor.cluster)) 
# corrplot(cor.cluster[hc$order, hc$order], order="original", method="square", type="full",tl.cex=0.3, tl.col="black")

# plot(hc)               # plot the dendrogram 
plot(hc, hang=-1)               # plot the dendrogram 

# cut tree into clusters
Ngroups = 6
rect.hclust(hc, k=Ngroups)
# get cluster IDs for the groups
groups <- cutree(hc, k=Ngroups)
groups.hc.order <- groups[hc$order]

dev.off()

# plot the groups
options$height = 3000
options$width = 3000

# Plot the clusters
for (k in 1:Ngroups){
  fname <- sprintf("%s_cluster_%s.png", method, k)
  png(filename=sprintf("../results/cluster/%s", fname), width=options$width, height=options$height, res=options$res)
  g <- groups.hc.order[groups.hc.order==k]
  N <- ceiling(sqrt(length(g)))
  par(mfrow=c(N,N))
  for (name in names(g)){
    f_single_plot(name_A=name) 
  }
  par(mfrow=c(1,1))  
  dev.off()
}