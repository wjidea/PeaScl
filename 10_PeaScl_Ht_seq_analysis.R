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
Interact <- c("L12","L12","L24","L24","L48","L48","P12","P12","P24","P24","P48","P48","S12","S12",
              "S24","S24", "S48","S48" )
## DESeq data set needs to be a data.frame
# SampleFile will connect independent files with the different conditions
samplename <- unlist(lapply(sampleFiles,function(x) gsub("\\.htout\\.new","",x)))


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
# colData(ddsHTSeq)$Hosts <- factor(colData(ddsHTSeq)$Hosts, levels=c("Media_CK","Lifter","PI240515"))
# colData(ddsHTSeq)$TP <- factor(colData(ddsHTSeq)$TP, levels=c("12","24","48"))
# colData(ddsHTSeq)$HostsTP <- factor(colData(ddsHTSeq)$, levels=c("12","24","48"))
# Normalize reads among gene expression
dds<-DESeq(ddsHTSeq)
dds$nested <- factor(c("L12","L12","L24","L24","L48","L48","P12","P12","P24","P24","P48","P48","S12","S12",
                       "S24","S24", "S48","S48"))
resall <- results(dds)
resall <- res[order(resall$padj),]
res <- results(dds,contrast = c("Hosts","PI240515","Lifter"), )
res <- res[order(res$padj),]
head(res)
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
plot(logcounts[,1], logcounts[,2], cex=.1)

# rlog transforamtion of data
# this takes ~15 seconds
rld <- rlogTransformation(dds)

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
plot(assay(rld)[,1], assay(rld)[,2], cex=.1)


resBigFC <- results(dds, cooksCutoff = 1, altHypothesis="greaterAbs")
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)

resSort <- res[order(res$pvalue),]
head(resSort)
k <- counts(dds)[rownames(resSort)[1],]
cond <- with(colData(dds), factor(paste(TP, Hosts)))
par(mar=c(8,5,2,2))
stripchart(log2(k + 1) ~ cond, method="jitter", vertical=TRUE, las=2)


res_time <- results(dds, contrast=c("TP","12","24"))
res_timeSort <- res_time[order(res_time$pvalue),]
head(res_timeSort,10)

res_Lifter <- results(dds, contrast=list("HostsLifter.TP12","HostsLifter.TP24"), pAdjustMethod = "BH")
res_timeSort <- res_time[order(res_time$padj),]
head(res_timeSort,10)



# #biocLite("org.Sc.sgd.db")
# library(org.Sc.sgd.db)
# keytypes(org.Sc.sgd.db)
# keytypes(org.Mm.eg.db)
# head(rownames(dds),20)
# geneinfo <- select(org.Sc.sgd.db, keys=rownames(resSort)[1:20],
#                    columns=c("PFAM","GO","EVIDENCE"), 
#                    keytype="PFAM")
# geneinfo
# 

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
