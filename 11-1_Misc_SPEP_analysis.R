# Install packages from Bioconductor
#source("http://bioconductor.org/biocLite.R")
# Dependencies
library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)
library(genefilter)
library(RColorBrewer)
# Function to draw plot based on gene name and contrast group
draw_genename <- function (resfunc,name){
  topGene <- grep(pattern = name,x = rownames(resfunc),value=TRUE)
  chartTitle <- paste("Gene Name: ",name)
  gdata <- plotCounts(dds , gene=topGene, intgroup=c("Hosts","TP"), 
                      returnData=TRUE,normalized = TRUE)
  plot_gene <- ggplot(gdata, aes(x=TP, y=count,color=Hosts, group=Hosts)) + 
    geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
  plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
}

draw_genename(assay(rld),"SS1G_12776")
test <- as.character(spep[spep$Gene_ID %in% "SS1G_12776",]$Annotation)


spep <- read.csv("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/SPEP.csv",header = TRUE,)
counts_normal <- as.data.frame(counts(dds,normalized=TRUE))
firstrow <- as.list(rownames(counts_normal))

GeneID_all <- unlist(lapply(firstrow,function(x) {
  ifelse(as.vector(unlist(strsplit(x,"\\|")))[2] %in% "Unknown",
         as.vector(unlist(strsplit(x,"\\|")))[1],
         as.vector(unlist(strsplit(x,"\\|")))[2])}))
counts_normal <- cbind(GeneID_all,counts_normal )

spepcount <- data.frame()
for (prot in spep$Gene_ID) {
  if (prot %in% GeneID_all) spepcount <- rbind(spepcount,counts_normal[counts_normal$GeneID_all %in% prot,])
}
spepcountSort <- spepcount[order(spepcount$GeneID_all),]
spepcountSubset <- subset(spepcountSort, select = -c(GeneID_all) )

# pc <- prcomp( t( spepcountSubset ) )
# plot(pc$x[,1], pc$x[,2], 
#      col=colData(dds)$TP, 
#      pch=as.numeric(colData(dds)$Hosts)+14)
# legend("bottomright",legend = c("Media-CK","Lifter","PI"),pch=c(15,16,17))
# legend("topright",legend=c("12","24","48"),col=c("Green","Red","Black"),pch=16)

# df here is the DE genes that belong to secreted proteins

df <- data.frame()
for (i in rownames(T24PI)){
  x<-as.vector(unlist(strsplit(i,"\\|")))[2]
  if (x %in% spepcountSort$GeneID_all){
    df <- rbind(df,spep[spep$Gene_ID %in% x,])
    df <- df[order(df$Effector_candidate),]
  }
}
draw_genename(T24PI,"SS1G_12776")
test <- as.character(spep[spep$Gene_ID %in% "SS1G_12776",]$Annotation)


a <- rownames(T12Lift)
b <- rownames(T24Lift)
c <- rownames(T12PI)
d <- rownames(T24PI)
res_FDR <- list(PLTP12,PLTP24,PLTP48,PICK12,PICK24,PICK48,LTCK12,
                LTCK24,LTCK48,T12Lift,T24Lift,T12PI,T24PI)

all_DENames <- lapply(res_FDR,function(x) rownames(x))
DENames_union <- Reduce(union,all_DENames)
write.table(DENames_union,sep="\t",file = "all_DENames.txt",quote = FALSE,row.names = FALSE)

# setA<-c("a", "b", "c", "d", "e")
# setB<-c("d", "e", "f", "g")
# setC<-c("f", "g","h","i")
# Reduce(union,list(setC,setA,setB))

##### Heat map for candidate effector genes ##########

topVarGenes <- head(order(-rowVars(assay(rld))),100)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
colors2 <- c("Blue","Yellow","Red")
sidecols <- c("grey","dodgerblue","lightgreen")[ rld$TP ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$Hosts,"-",rld$TP)
# heatmap.2(mat, trace="none", , ColSideColors=sidecols,
#           labRow=TRUE, key=FALSE,lhei=c(2,8), lwid=c(1,4))
heatmap(mat,col=colors,cexRow = 0.5,ColSideColors = sidecols,scale="row",)
dev.off()

spep_effector <- spep[spep$Effector_candidate %in% c("YES"),]
spep_effector_count <- data.frame()
for (prot in spep_effector$Gene_ID) {
  if (prot %in% counts_normal$GeneID_all) spep_effector_count <- rbind(spep_effector_count,counts_normal[counts_normal$GeneID_all %in% prot,])
}
spep_effector_count<- subset(spep_effector_count, select = -c(GeneID_all))
spep_effector_count <- as.matrix(spep_effector_count)
spep_effector_count <- spep_effector_count - rowMeans(spep_effector_count)
colnames(spep_effector_count) <- paste0(rld$Hosts,"-",rld$TP)
rownames(spep_effector_count)
top100 <- head(order(-rowVars(spep_effector_count)),20)
top_effector_count <- spep_effector_count[top100,]
par(mar=c(0.1,0.1,0.1,2))
heatmap.2(top_effector_count, trace="none", dendrogram="col", col=colors, ColSideColors=sidecols,
          cexCol=1, lhei=c(1.5,2), lwid=c(1,6), key=TRUE, keysize=0.02,revC=FALSE)
heatmap(top_effector_count,col=colors,scale = "none")
