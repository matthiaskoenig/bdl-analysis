# Gene set enrichment analysis


Overrepresentation analysis within the clusters.
One of the main uses of the GO is to perform enrichment analysis on gene sets. For example, given a set of genes that are up-regulated under certain conditions, an enrichment analysis will find which GO terms are over-represented (or under-represented) using annotations for that gene set.

Three subsets of GO, biological processes, cellular components and molecular function.
``` {r}
uniprot_from_name <- function(name){
  return(ProbeInformation(name)$Entry)
}


hc.res <- hclust.res[[method]]  # clustering for correlation method
groups.hc.order <- hc.res$groups.hc.order
for (k in 1:Ngroups){
  g <- groups.hc.order[groups.hc.order==k]
  print('---------------------------------')
  print(sprintf("Cluster %s", k))
  print(names(g))  # names for the cluster representatives
  uniprot <- lapply(names(g), uniprot_from_name)
  names(uniprot) <- names(g)
  # filter the non-existing UniProts, i.e. for non-gene data
  uniprot <- Filter(Negate(is.null), uniprot)
  print(as.data.frame(uniprot))
}

install.packages("clusterProfiler")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("topGO")

# make the topGO object
levels(BDLfactors$ftype)

all.genes <- BDLfactors[BDLfactors$ftype %in% c("GE_ADME", "GE_Cytokines", "GE_Fibrosis"),]$id
all.genes[all.genes]
names(all.genes) <- all.genes

names(g)

library("topGO")
GOdata.BP <- new("topGOdata", ontology='BP', 
                 allGenes = all.genes, 
                 geneSel = names(g),
                 nodeSize = 1,
                 annotationFun = annFUN.db, affyLib=affyLib)


# set of all genes, i.e. the genes on all fluidigm chips
BDLfactors$ftype


```