library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)
# Prepare dataset
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output")
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



##### WGCNA Analysis ####
library(WGCNA)
library(matrixStats)
library(RColorBrewer)
#Important command, DO NOT OMIT
options(stringsAsFactors = FALSE)

#Reading data
data_expr0 <- as.data.frame(assay(rld))
#Sets for Lifter, PI and Media
nSets = 3
#WGCNA requires gene names on columns
data_expr0 <- t(data_expr0)

#Break data into sets
Lifter_expr <- data_expr0[1:6,]
PI_expr <- data_expr0[7:12,]
SclCK_expr <- data_expr0[13:18,]

#Remove columns that sum 0
Lifter_expr <- Lifter_expr[,colSums(Lifter_expr) != 0]
PI_expr <- PI_expr[,colSums(PI_expr) != 0]
SclCK_expr <- SclCK_expr[,colSums(SclCK_expr) != 0]

#Calculate coefficient of variation
# Lifter_CV <- colMeans(Lifter_expr)/colSds(Lifter_expr)
# count(Lifter_CV, value = 'NA')
# 
# PI_CV <- colSds(PI_expr,na.rm=TRUE)/colMeans(PI_expr,na.rm = TRUE)
# count(PI_CV, value =NA)
# 
# SclCK_CV <- colSds(SclCK_expr,na.rm=TRUE)/colMeans(SclCK_expr,na.rm = TRUE)
# count(SclCK_CV, value = 'NA')
 

#Set data
setLabels <- c("Lifter-Scl","PI240515-Scl","Scl-CK")
shortLabels <- c("Lifter","PI","Scl")
#Form multiset expression data
multiExpr <- vector(mode="list", length=nSets)

multiExpr[[1]] <- list(data = Lifter_expr)
multiExpr[[2]] <- list(data = PI_expr)
multiExpr[[3]] <- list(data = SclCK_expr)
#CheckSets for samples and structure nGenes=14318
checkSets(multiExpr)

#Checking missing values, if $allOK data can be used.
gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK

#Remove outliers and cleaning
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
  printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
  for (set in 1:nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
}

#Check sets again nGenes=12796
checkSets(multiExpr)

##Clustering samples
#Loop to generate clustering
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

#Plot trees for clustering
#No outliers present based on clustering
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for(set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2)[[2]])
  }
collectGarbage()

# Plot the results:
colors=brewer.pal(5,"Set1")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
#sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
}

### Calculation of network adjacencies
softPower = 6
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower

### Calculation of topological overlap
# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);

