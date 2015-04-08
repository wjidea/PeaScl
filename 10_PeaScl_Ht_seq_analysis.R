library(DESeq2)
library(Biobase)
library(GenomicRanges)
#Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_scl/Broad_gff/")
directory <- "/home/peascl/scl_de/htout_STAR_scl/Broad_gff/"
sampleFiles <- grep("htout",list.files(directory),value=TRUE)
HostNames <- c(rep("Lifter",6),rep("PI240515",6),rep("Media_CK",6))
TimePoint <- c("12","12","24","24","48","48",
               "12","12","24","24","48","48",
               "12","12","24","24","48","48")
#DESeq data set needs to be a data.frame
sampleTable <- data.frame(sampleName=sampleFiles,
                          fileName=sampleFiles,
                          Hosts=HostNames,
                          TP=TimePoint
)
#Read in dataset
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, 
                                     directory=directory, 
                                     design=~Hosts + TP + Hosts:TP)
# relevel to make sure the order of the parameters in the right order
# otherwise relvel is needed
colData(ddsHTSeq)$Hosts <- factor(colData(ddsHTSeq)$Hosts, levels=c("Media_CK","Lifter","PI240515"))
colData(ddsHTSeq)$TP <- factor(colData(ddsHTSeq)$TP, levels=c("12","24","48"))
# normalize reads among gene expression
dds<-DESeq(ddsHTSeq)
resall <- results(dds)
resall <- res[order(resall$padj),]
res <- results(dds,contrast = c("Hosts","PI240515","Lifter"), )
res <- res[order(res$padj),]
head(res)
plotMA(dds,ylim=c(-2,2))

dds <- estimateSizeFactors( dds )
sizeFactors(dds)
colSums(counts(dds))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )

plot(pc$x[,1], pc$x[,2], 
     col=colData(dds)$TP, 
     pch=as.numeric(colData(dds)$Hosts)+14)
legend("bottomright",legend = c("Media-CK","Lifter","PI"),pch=c(15,16,17))
legend("topright",legend=c("12","24","48"),col=c("Green","Red","Black"),pch=16)
plot(hclust(dist(t(logcounts))), labels=colData(dds)$Hosts)
plot(hclust(dist(t(logcounts))), labels=colData(dds)$TP)
plot(logcounts[,1], logcounts[,2], cex=.1)

# this takes ~15 seconds
rld <- rlogTransformation(dds)
pc2 <- prcomp( t( assay(rld) ) )

plot(pc2$x[,1], pc2$x[,2],
     col=colData(rld)$TP, 
     pch=as.numeric(colData(rld)$Hosts)+14)
legend("topright",legend = c("Media-CK","Lifter","PI"),pch=c(15,16,17))
legend("bottomright",legend=c("12","24","48"),col=c("Green","Red","Black"),pch=16)
par(mfrow=c(1,2))
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$Hosts)
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$TP)
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


biocLite("org.Sc.sgd.db")
library(org.Sc.sgd.db)
keytypes(org.Sc.sgd.db)
keytypes(org.Mm.eg.db)
head(rownames(dds),20)
geneinfo <- select(org.Sc.sgd.db, keys=rownames(resSort)[1:20],
                   columns=c("PFAM","GO","EVIDENCE"), 
                   keytype="PFAM")
geneinfo

