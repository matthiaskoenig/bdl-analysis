# Version info: R 2.14.1, Biobase 2.15.3, GEOquery 2.23.2, limma 3.10.1
# R scripts generated  Sun Aug 23 10:24:35 EDT 2015

################################################################
#   Differential expression analysis with limma
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("GEOquery", "limma")
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE29776", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL81", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
sampleNames(gset)

# group names for all samples
sml <- c("BDL","BDL","BDL","sham","sham","sham");

# Under the top genes should be: Serpine1, Mr1, Pai1, Planh1

# Normalisation
biocLite("simpleaffy")
biocLite("affyPLM")
library(simpleaffy)
library(affyPLM)
gset.loess <- normalize.ExpressionSet.loess(gset)
gset.normalized <- gset.loess



# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset.normalized)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(BDL-sham, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)


# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=12483)
# tT[c("ID", "Gene.Symbol", "ENTREZ_GENE_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# load NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
write.table(tT, file="Wang2005.txt", row.names=F, sep="\t")


head(tT, 300)

tT[tT$Gene.symbol=="Mr1" , ]

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE29776", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL81", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
sml <- c("G0","G0","G0","G1","G1","G1")

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("BDL","sham")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE29776", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
