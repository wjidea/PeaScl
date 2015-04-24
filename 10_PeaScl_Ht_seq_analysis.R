library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)
# Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF/")
directory <- "/home/peascl/scl_de/htout_STAR_mergeGTF"
sampleFiles <- grep("htout",list.files(directory),value=TRUE)
HostNames <- c(rep("Lifter",6),rep("PI240515",6),rep("Media_CK",6))
TimePoint <- c("12","12","24","24","48","48",
               "12","12","24","24","48","48",
               "12","12","24","24","48","48")
## DESeq data set needs to be a data.frame
# SampleFile will connect independent files with the different conditions
samplename <- unlist(lapply(sampleFiles,function(x) gsub("\\.htout","",x)))
# prepare sample table
sampleTable <- data.frame(sampleName=samplename,
                          fileName=sampleFiles,
                          Hosts=HostNames,
                          TP=TimePoint
)
## Read in dataset
# SampleTable has all the inforamtion for files and conditions
# Directory contains the actual files
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                     directory=directory, 
                                     design=~Hosts*TP)

# Re-level to make sure the order of the parameters in the right order
# Otherwise re-level is needed
colData(ddsHTSeq)$Hosts <- factor(colData(ddsHTSeq)$Hosts, levels=c("Media_CK","Lifter","PI240515"))
colData(ddsHTSeq)$TP <- factor(colData(ddsHTSeq)$TP, levels=c("12","24","48"))

# Normalize reads among gene expression
dds<-DESeq(ddsHTSeq,test = "Wald")

# explore the result in general, and generate FDR
#resall <- results(dds,pAdjustMethod = "BH")
# sort by FDR
#resall <- res[order(resall$padj),]
#plotMA(dds,ylim=c(-2,2))

# Estimate size factors and plot 
dds <- estimateSizeFactors( dds )
sizeFactors(dds)
colSums(counts(dds))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

### Principal component analysis of samples ##
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )

## Before rlog transformation (variance stabilizing transformation)
# PCA
plot(pc$x[,1], pc$x[,2], 
     col=colData(dds)$TP, 
     pch=as.numeric(colData(dds)$Hosts)+14)
legend("bottomright",legend = c("Media-CK","Lifter","PI"),pch=c(15,16,17))
legend("topright",legend=c("12","24","48"),col=c("Green","Red","Black"),pch=16)
# Hierarchical clustering
plot(hclust(dist(t(logcounts))), labels=colData(dds)$Hosts)
plot(hclust(dist(t(logcounts))), labels=colData(dds)$TP)
# logcounts
# plot(logcounts[,1], logcounts[,3], cex=.1)

# rlog transforamtion of data
# this takes ~15 seconds
rld <- rlogTransformation(dds,blind=FALSE)

## PCA after rlog transformation of data
pc2 <- prcomp( t( assay(rld) ) )

# PCA plot
plot(pc2$x[,1], pc2$x[,2],
     col=colData(rld)$TP, 
     pch=as.numeric(colData(rld)$Hosts)+14)
legend("topright",legend = c("Media-CK","Lifter","PI"),pch=c(15,16,17))
legend("bottomright",legend=c("12","24","48"),col=c("Green","Red","Black"),pch=16)

# Hierarchical clustering
par(mfrow=c(1,2))
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$Hosts)
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$TP)

# log counts
par(mfrow=c(1,1))
plot(assay(rld)[,1], assay(rld)[,3], cex=.1)

# show results with the genes that are up or down regulated
resBigFC <- results(dds, altHypothesis="greaterAbs",pAdjustMethod = "BH",)
resBigFCsort <- resBigFC[order(resBigFC$padj),]
head(resBigFCsort)
Gene <- rownames(resBigFCsort)[1]
plotMA(resBigFC, ylim=c(-3,3))
#abline(h=c(-1,1),lwd=5)

#calculate contrast with significant threshold with FDR significance level of .1
Cal.ContrastSig <- function(vec1, vec2,lV=c(1,-1),threshold=0.1){
  res_func <- results(dds, contrast = list(vec1,vec2),
                      pAdjustMethod = "BH",listValues=lV)
  res_funcSort <- res_func[order(res_func$padj),]
  resSig <- subset(res_funcSort, padj < threshold)
  resSig
}
Cal.ContrastSig(c("TP12","HostsLifter.TP12"),c("TP24","HostsLifter.TP24"),c(1,-1),0.1)

library("RColorBrewer")
library("gplots")
select <- order(rowVars(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)


heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = TRUE, scale="row",
          trace="none", margin=c(10, 2), 
          lhei=c(2,5), lwid=c(1,5), keysize=0.1)

distsRL <- dist(t(assay(rld))) 
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=" : "))
#if you just want the conditions use this line
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16),symm = FALSE)

print(plotPCA(rld,intgroup =c("TP")))

res_lift_CK <- results(dds,contrast = list("HostsLifter.TP12", "HostsMedia_CK.TP12"))
res_lift_CKSort <- res_lift_CK[order(res_lift_CK$padj),]
head(res_lift_CKSort,50)

res_lift_CK <- results(dds,contrast = list("HostsMedia_CK.TP12","HostsLifter.TP12"))
res_lift_CKSort <- res_lift_CK[order(res_lift_CK$padj),]
head(res_lift_CKSort,10)