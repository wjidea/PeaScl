library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)

# Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/")
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
                          TP=TimePoint)
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

# Estimate size factors and plot 
dds <- estimateSizeFactors( dds )

# rlog transforamtion of data
rld <- rlogTransformation(dds,blind=FALSE,)

########### End of data DEseq2 prep #################

# PCA after rlog transformation of data
pc2 <- prcomp( t( assay(rld) ) )
plotPCA(assay(rld),intgroup = c(Hosts,TP))
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
with(resBigFCsort[Gene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, Gene, 
       pos=4, col="dodgerblue")
})

####################### TIME Effect  ###########################
# Detect time effect in host Lifter and PI
T12Lift <- Cal.ContrastSig(c("TP12","HostsLifter.TP12"),c("TP24","HostsLifter.TP24"))
T24Lift <- Cal.ContrastSig(c("TP24","HostsLifter.TP24"),c("TP48","HostsLifter.TP48"))
T12PI <- Cal.ContrastSig(c("TP12","HostsPI240515.TP12"),c("TP24","HostsPI240515.TP24"))
T24PI <- Cal.ContrastSig(c("TP24","HostsPI240515.TP24"),c("TP48","HostsPI240515.TP48"))
Time_effects <- list(T12Lift,T24Lift,T12PI,T24PI)
#lapply(Time_effects,function(x) write.table(x,"Time_Effect",append=TRUE,quote=FALSE,sep = ",",col.names =TRUE ) )
####################### HOST Effect  ###########################
# Lifter PI HOST effect in TP12
PLTP12 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP12"),c("HostsLifter","HostsLifter.TP12"))
PLTP24 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP24"),c("HostsLifter","HostsLifter.TP24"))
PLTP48 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP48"),c("HostsLifter","HostsLifter.TP48"))
Hosts_effects <- list(PLTP12,PLTP24,PLTP48)
#lapply(Hosts_effects,function(x) write.table(x,"Hosts_effects",append=TRUE,quote=FALSE,sep = ",",col.names =TRUE ) )

####################### PI Compare to Media #######################
PICK12 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP12"),c("HostsMedia_CK","HostsMedia_CK.TP12"))
PICK24 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP24"),c("HostsMedia_CK","HostsMedia_CK.TP24"))
PICK48 <- Cal.ContrastSig(c("HostsPI240515","HostsPI240515.TP48"),c("HostsMedia_CK","HostsMedia_CK.TP48"))
PIMedia_effects <- list(PLTP12,PLTP24,PLTP48)
#lapply(PIMedia_effects,function(x) write.table(x,"PIMedia_effects",append=TRUE,quote=FALSE,sep = ",",col.names =TRUE ) )

####################### Lifter Compare to Media #######################
LTCK12 <- Cal.ContrastSig(c("HostsLifter","HostsLifter.TP12"),c("HostsMedia_CK","HostsMedia_CK.TP12"))
LTCK24 <- Cal.ContrastSig(c("HostsLifter","HostsLifter.TP24"),c("HostsMedia_CK","HostsMedia_CK.TP24"))
LTCK48 <- Cal.ContrastSig(c("HostsLifter","HostsLifter.TP48"),c("HostsMedia_CK","HostsMedia_CK.TP48"))
LTMedia_effect <- list(LTCK12,LTCK24,LTCK48)
#lapply(LTMedia_effect,function(x) write.table(x,"LTMedia_effect",append=TRUE,quote=FALSE,sep = ",",col.names =TRUE ) )

####################### All Compare to Media #######################
Trt_CK12 <- Cal.ContrastSig(c("HostsPI240515.TP12","HostsLifter.TP12"),c("HostsMedia_CK.TP12"),lV=c(1/2,-1))
Trt_CK24 <- Cal.ContrastSig(c("HostsPI240515.TP24","HostsLifter.TP24"),c("HostsMedia_CK.TP24"),lV=c(1/2,-1))
Trt_CK48 <- Cal.ContrastSig(c("HostsPI240515.TP48","HostsLifter.TP48"),c("HostsMedia_CK.TP48"),lV=c(1/2,-1))
PlantVsMedia <- list(Trt_CK12,Trt_CK24,Trt_CK48)
#lapply(PlantVsMedia,function(x) write.table(x,"PlantVsMedia",append=TRUE,quote=FALSE,sep = ",",col.names =TRUE ) )


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