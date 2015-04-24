# Install 
#install.packages("gProfileR")

# Dependencies
library(DESeq2)
library(Biobase)
library(GenomicRanges)
library(ggplot2)
library(gplots)
library(genefilter)
library(RColorBrewer)
library(biomaRt)
library(gProfileR)

#calculate contrast with significant threshold with FDR significance level of .1
Cal.ContrastSig <- function(vec1, vec2,lV=c(1,-1),threshold=0.1){
  res_func <- results(dds, contrast = list(vec1,vec2),
                      pAdjustMethod = "BH",listValues=lV)
  res_funcSort <- res_func[order(res_func$padj),]
  resSig <- subset(res_funcSort, padj < threshold)
  resSig
}

# Function to draw plot based on gene name and contrast group
draw_genename <- function (gene,title){
  chartTitle <- paste("Gene Name: ",title)
  gdata <- plotCounts(dds,gene=gene, intgroup=c("Hosts","TP"), 
                      returnData=TRUE,normalized = TRUE)
  plot_gene <- ggplot(gdata, aes(x=TP, y=count,color=Hosts, group=Hosts)) + 
    geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
  plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
}

# Function to draw plot based on gene name with rld transformed data
draw_genename_rld <- function (gene,title){
  chartTitle <- paste("Gene Name: ",title)
  data <- assay(rld)[rownames(assay(rld))==gene,]
  single_gene_rld <- as.data.frame(data)
  single_gene_rld$Hosts <- paste0(rld$Hosts)
  single_gene_rld$TP <- paste0(rld$TP)
  plot_gene <- ggplot(single_gene_rld, aes(x=TP, y=data,color=Hosts, group=Hosts)) + 
    geom_point(position=position_jitter(width=0.1,height=0),size=3,shape=16) 
  plot_gene + labs(title = chartTitle) + stat_smooth(se=F,method="loess")
}

# Function to draw pie chart of the Go enrichment result from cutree function
draw_pie <- function(clade,ti="clade"){
  pie_v <- clade$V3
  lbls <- paste(clade$V2,"\n",clade$V5)
  main_title <- paste(ti,"protein with annotations:",sum(clade$V3))
  pie(pie_v,col=brewer.pal(length(pie_v), "Set3"),labels = lbls,main = main_title)
}
