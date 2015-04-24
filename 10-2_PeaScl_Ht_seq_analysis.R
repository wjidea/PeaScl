#Working directory
setwd("/home/peascl/scl_de/htout_STAR_mergeGTF/")

##Source previous data
source("/home/peascl/scl_de/Exec/10_PeaScl_Ht_seq_analysis.R")

##### WGCNA Analysis ####
library(WGCNA)
library(matrixStats)
library(RColorBrewer)
library(car)

#Important command, DO NOT OMIT
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
Sys.setenv(OMP_NUM_THREADS=4)

#Reading data
data_expr0 <- as.data.frame(assay(rld))

#Filtering data
histinfo <- hist(rowVars(as.matrix(data_expr0)),breaks = 100)
rug(rowVars(as.matrix(data_expr0)))
data_expr0 <- data_expr0[rowVars(as.matrix(data_expr0))>=0.5,]

#WGCNA requires gene names on columns
#number of genes 4767
data_expr0 <- t(data_expr0)  
  
#Checking missing values, if $allOK data can be used.
gsg <- goodSamplesGenes(data_expr0, verbose = 3)
gsg$allOK

#Remove outliers and cleaning
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data_expr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data_expr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data_expr0 = data_expr0[gsg$goodSamples, gsg$goodGenes]
}

#Check number of genes again Genes=4767
dim(data_expr0)


##Clustering samples
#hclust of samples
sampleTree <- hclust(dist(data_expr0), method = "average")

#Plot trees for clustering
#No outliers present based on clustering
plot(sampleTree, main = "Sample clustering to detect outliers",
       xlab="", sub="", cex.lab = 1, cex.axis = 1, cex.main = 1.5)
#saving data for final analysis
data_expr <- data_expr0

###Associate sample data to expression values
#Data for samples
sample_data <- sampleTable
#Recoding data media=10, lifter=20 and PI=30
sample_data$Hosts <- recode(sample_data$Hosts, "'Media_CK'=10;'Lifter'=20;'PI240515'=30")
sample_data$Host_TP <- paste(sample_data$Hosts, sample_data$TP, sep = "")
sample_data <- sample_data[,c(1,3,4,5)]


#Match samples to sample data
samples <- row.names(data_expr)
dat_row <- match(samples, sample_data$sampleName)
datSample <- sample_data[dat_row, -1]
row.names(datSample) <- sample_data[dat_row,1]

datSample$TP <- as.numeric(datSample$TP)
datSample$Host_TP <- as.numeric(datSample$Host_TP)

collectGarbage()

###re-cluster sampel to correlate sample data
sampleTree2 <- hclust(dist(data_expr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datSample, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datSample),
                    main = "Sample dendrogram and trait heatmap")

#####Step by step network construction
### 1. Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))

#Call topology analysis function
sft <- pickSoftThreshold(data_expr, powerVector=powers, verbose = 5)

# plot the returned scale free analysis tables
#sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### 2. Calculation of network similarity and adjacencies
softPower = 24
adjacency <- adjacency(data_expr, power = softPower)

### 3. Calculation of topological overlap matrix
tom <- TOMsimilarity(adjacency)
dissTOM <- 1 - tom

### 4. Clustering using TOM
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
#This creates a total of 55 different modules
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

### 5. Merge gene modules with similar expression
# Calculate eigengenes
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.20
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(data_expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs

#New clusters
#sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "PeaScl_WGCNA-stepByStep.RData")

##### Associating modules to sample data
### Quatifying module-trait associations
# Define numbers of genes and samples
nGenes <- ncol(data_expr)
nSamples <- nrow(data_expr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(data_expr, moduleColors)$eigengenes
#Sort eigenvectors by similarity
MEs2 <- orderMEs(MEs0)
#Correlations using Pearson
moduleTraitCor <- cor(MEs2, datSample, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 9, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datSample),
               yLabels = names(MEs2),
               ySymbols = names(MEs2),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## Gene relationship and important modules
# Define variable weight containing the H_TP column of datSample
H_TP <- as.data.frame(datSample$Host_TP)
names(H_TP) = "host_TP"
# names (colors) of the modules
modNames <- substring(names(MEs2), 3)
geneModuleMembership <- as.data.frame(cor(data_expr, MEs2, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance <- as.data.frame(cor(data_expr, H_TP, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(H_TP), sep="");
names(GSPvalue) = paste("p.GS.", names(H_TP), sep="");

### GS and MM
module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Host",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

####Visualize the gene network
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM <- dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA
# Call the plot function
#sizeGrWindow(9,9)
######Error here: Error: C stack usage  7971408 is too close to the limit######
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
######

####Visualize network of eigengenes
# Recalculate module eigengenes
MEs3 <- moduleEigengenes(data_expr, moduleColors)$eigengenes
# Isolate Hosts from the sample data
Hosts <- as.data.frame(datSample$Hosts)
names(Hosts) <- "Hosts"
TP <- as.data.frame(datSample$TP)
names(TP) <- "Time"
HTP <- as.data.frame(datSample$Host_TP)
names(HTP) <- "Host-TP"
# Add the Hosts to existing module eigengenes
MET <- orderMEs(cbind(MEs3, Hosts,TP,HTP))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle= 90)

#### Exporting to cytoscape
# Select modules
modules = c("lightcyan", "saddlebrown","grey60","cyan","magenta");
# Select module probes
probes <- colnames(data_expr)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]
modGenes <- annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM <- tom[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

#Extract gene names
write.table(colnames(data_expr)[moduleColors=="lightcyan"], file="lightcyan_DR.txt", quote = FALSE,row.names = FALSE)
write.table(colnames(data_expr)[moduleColors=="saddlebrown"], file="saddlebrown_DR.txt", quote = FALSE,row.names = FALSE)
write.table(colnames(data_expr)[moduleColors=="grey60"], file="grey60_DR.txt", quote = FALSE,row.names = FALSE)
write.table(colnames(data_expr)[moduleColors=="cyan"], file="cyan_UR.txt", quote = FALSE,row.names = FALSE)
write.table(colnames(data_expr)[moduleColors=="magenta"], file="magenta_UR.txt", quote = FALSE,row.names = FALSE)
#extract gene_ID
write.table(modProbes,file = "Gene_ID.txt",quote = FALSE,row.names = FALSE)

# test1 <- c("SS1G_00012")
# 
# draw_genename_rld <- function (gene,title){
#   chartTitle <- paste("Gene Name: ",title)
#   data <- assay(rld)[rownames(assay(rld))==gene,]
#   single_gene_rld <- as.data.frame(data)
#   single_gene_rld$Hosts <- paste0(rld$Hosts)
#   single_gene_rld$TP <- paste0(rld$TP)
#   plot_gene <- ggplot(single_gene_rld, aes(x=TP, y=data,color=Hosts, group=Hosts)) + 
#     geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
#   plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
# }