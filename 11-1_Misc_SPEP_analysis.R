# Dependencies
library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)

# Function to draw plot based on gene name and contrast group
draw_genename <- function (resfunc,name){
  topGene <- grep(pattern = name,x = rownames(resfunc),value=TRUE)
  chartTitle <- paste("Gene Name: ",name)
  gdata <- plotCounts(dds , gene=topGene, intgroup=c("Hosts","TP"), 
                      returnData=TRUE,normalized = TRUE)
  plot_gene <- ggplot(gdata, aes(x=TP, y=count,color=Hosts, group=Hosts)) + 
    geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
  p <- plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
  return(list(p,as.vector(unlist(strsplit(topGene,"\\|")))[2]))
}




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
for (i in rownames(T24Lift)){
  x<-as.vector(unlist(strsplit(i,"\\|")))[2]
  if (x %in% spepcountSort$GeneID_all){
    df <- rbind(df,spep[spep$Gene_ID %in% x,])
    df <- df[order(df$Effector_candidate),]
  }
}
draw_genename(T24Lift,"SS1G_12024")


