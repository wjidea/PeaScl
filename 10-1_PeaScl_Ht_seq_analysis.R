library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)
# Install necessary packages
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF")
directory <- "/home/peascl/scl_de/htout_STAR_mergeGTF"
sampleFiles <- grep("htout",list.files(directory),value=TRUE)
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
Gene <- rownames(resBigFCsort)[1]
plotMA(resBigFC, ylim=c(-3,3))
#abline(h=c(-1,1),lwd=5)
with(resBigFCsort[Gene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, as.vector(unlist(strsplit(Gene,"\\|")))[2], 
       pos=4, col="dodgerblue")
})

# Cal.contrast <- function(vec1, vec2,lV=c(1,-1),top=50){
#   res_func <- results(dds, contrast = list(vec1,vec2),
#                       pAdjustMethod = "BH",listValues=lV)
#   res_funcSort <- res_func[order(res_func$padj),]
#   head(res_funcSort,top)
# }
  # Define conditions
# cond <- with(colData(dds), factor(paste(Hosts,TP)))
#   # VERY good function to show gene expression level between samples 
#   # Give name to each column in the stripchart
#   # obtain normalized read from the top DE genes 
#   # [2] is the rank of the gene in the dds sorted result (i.e. resSort)
# draw_stripchart <- function(x,n){
#   k <- counts(dds,normalized=TRUE)[rownames(x)[n],]
#   par(mar=c(8,5,2,2))
#   stripchart(log2(k + 1) ~ cond, method="jitter", vertical=TRUE, las=2,)
#   print(rownames(x)[n])
# }
#calculate contrast with significant threshold with FDR significance level of .1
Cal.ContrastSig <- function(vec1, vec2,lV=c(1,-1),threshold=0.1){
  res_func <- results(dds, contrast = list(vec1,vec2),
                      pAdjustMethod = "BH",listValues=lV)
  res_funcSort <- res_func[order(res_func$padj),]
  resSig <- subset(res_funcSort, padj < threshold)
  resSig
}
#Cal.ContrastSig(c("TP12","HostsLifter.TP12"),c("TP24","HostsLifter.TP24"),c(1,-1),0.1)

draw_geneplot <- function (resfunc,n=1){
  topGene <- rownames(resfunc)[n]
  chartTitle <- paste("Gene Name: ",as.vector(unlist(strsplit(topGene,"\\|")))[2])
  gdata <- plotCounts(dds , gene=topGene, intgroup=c("Hosts","TP"), 
                      returnData=TRUE,normalized = TRUE)
  plot_gene <- ggplot(gdata, aes(x=TP, y=count,color=Hosts, group=Hosts)) + 
    geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
  p <- plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
  return(list(p,as.vector(unlist(strsplit(topGene,"\\|")))[2]))
}

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

# output summary vector



#draw_stripchart(Trt_CK,10)
res_test <- results(dds, contrast = list(c("TP12","HostsLifter.TP12"),c("TP24","HostsLifter.TP24")),pAdjustMethod = "BH",listValues=c(1,-1))
summary(Trt_CK24)

sum(res_test$padj < 0.1, na.rm=TRUE)
draw_geneplot(PICK12,1)


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