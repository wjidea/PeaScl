library(DESeq2)
library(Biobase)
library(GenomicRanges)
# Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF")
directory <- "/home/peascl/scl_de/htout_STAR_mergeGTF"
sampleFiles <- grep("htout.new",list.files(directory),value=TRUE)
HostNames <- c(rep("Lifter",6),rep("PI240515",6),rep("Media_CK",6))
TimePoint <- c("12","12","24","24","48","48",
               "12","12","24","24","48","48",
               "12","12","24","24","48","48")
## DESeq data set needs to be a data.frame
# SampleFile will connect independent files with the different conditions
samplename <- unlist(lapply(sampleFiles,function(x) gsub("\\.htout\\.new","",x)))
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

#if relevel is needed
#dds$Hosts <- relevel(dds$Hosts, "Media_CK")

# explore the result in general, and generate FDR
resall <- results(dds,pAdjustMethod = "BH")
# sort by FDR
resall <- res[order(resall$padj),]
plotMA(dds,ylim=c(-2,2))

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
rld <- rlogTransformation(dds,blind=FALSE,)

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
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)


# generate results with contrast
# res.contrast <- function =(var1,var2,var3,var4){
# }

####################### TIME Effect  ###########################
# 12 - 24 time effect in host Lifter
res_1 <- results(dds,contrast = list(
                        c("TP12","HostsLifter.TP12"), 
                        c("TP24","HostsLifter.TP24")),
                        pAdjustMethod = "BH")
res_1_Sort <- res_1[order(res_1$padj),]
T12Lift <- head(res_1_Sort,50)
Test_res <- Cal.contrast("HostsPI240515","HostsPI240515.TP24","HostsLifter","HostsLifter.TP24")
# 24 - 48 time effect in host Lifter
res_2 <- results(dds,contrast = list(
  c("TP24","HostsLifter.TP24"), 
  c("TP48","HostsLifter.TP48")),
  pAdjustMethod = "BH")
res_2_Sort <- res_2[order(res_2$padj),]
T24Lift <- head(res_2_Sort,50)

# 12 - 24 time effect in host PI
res_3 <- results(dds,contrast = list(
  c("TP12","HostsPI240515.TP12"), 
  c("TP24","HostsPI240515.TP24")),
  pAdjustMethod = "BH")
res_3_Sort <- res_3[order(res_3$padj),]
T12PI <- head(res_3_Sort,50)

# 24 - 48 time effect in host PI
res_4 <- results(dds,contrast = list(
  c("TP24","HostsPI240515.TP24"), 
  c("TP48","HostsPI240515.TP48")),
  pAdjustMethod = "BH")
res_4_Sort <- res_4[order(res_4$padj),]
T24PI <- head(res_4_Sort,50)


####################### HOST Effect  ###########################
# Lifter PI HOST effect in TP12
res_5 <- results(dds,contrast = list(
  c("HostsPI240515","HostsPI240515.TP12"), 
  c("HostsLifter","HostsLifter.TP12")),
  pAdjustMethod = "BH")
res_5_Sort <- res_5[order(res_5$padj),]
PLTP12 <- head(res_5_Sort,50)

# Lifter PI HOST effect in TP24
res_6 <- results(dds,contrast = list(
  c("HostsPI240515","HostsPI240515.TP24"), 
  c("HostsLifter","HostsLifter.TP24")),
  pAdjustMethod = "BH")
res_6_Sort <- res_6[order(res_6$padj),]
PLTP24 <- head(res_6_Sort,50)

Cal.contrast <- function(var1, var2, var3, var4){
  res_func <- results(dds, contrast = list(
    c(var1,var2),c(var3,var4)
  ),pAdjustMethod = "BH")
  res_funcSort <- res_func[order(res_func$padj),]
  head(res_funcSort,50)
}

Test_res <- Cal.contrast("HostsPI240515","HostsPI240515.TP24","HostsLifter","HostsLifter.TP24")




draw_stripchart(res_1_Sort,3)

cond <- with(colData(dds), factor(paste(Hosts,TP)))
draw_stripchart <- function(x,n){
# VERY good function to show gene expression level between samples 
# Give name to each column in the stripchart
# obtain normalized read from the top DE genes 
# [2] is the rank of the gene in the dds sorted result (i.e. resSort)
k <- counts(dds,normalized=TRUE)[rownames(x)[n],]
par(mar=c(8,5,2,2))
stripchart(log2(k + 1) ~ cond, method="jitter", vertical=TRUE, las=2,)
}

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

distsRL <- dist(t(assay(rld))) 
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=" : "))
#if you just want the conditions use this line
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16),symm = FALSE)

print(plotPCA(rld,intgroup =c("Hosts","TP") ))

res_lift_CK <- results(dds,contrast = list("HostsLifter.TP12", "HostsMedia_CK.TP12"))
res_lift_CKSort <- res_lift_CK[order(res_lift_CK$padj),]
head(res_lift_CKSort,50)

res_lift_CK <- results(dds,contrast = list("HostsMedia_CK.TP12","HostsLifter.TP12"))
res_lift_CKSort <- res_lift_CK[order(res_lift_CK$padj),]
head(res_lift_CKSort,10)