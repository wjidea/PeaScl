source("/home/peascl/scl_de/Exec/Functions_dependencies.R")

# Code to draw genes of interest temporal expression plots
draw_genename(gene = "SS1G_05875",title = "Pac1")
draw_genename(gene = "SS1G_08814",title = "Ss-Odc1")
draw_genename(gene = "SS1G_10796",title = "Ss-Odc2")
draw_genename(gene = "SS1G_00004",title = "GOI")
draw_genename(gene = "SS1G_11992",title = "GOI")
draw_genename_rld(gene = "SS1G_11992",title = "GOI")

draw_genename_rld(gene = "SS1G_02553",title = "SS1G_02553")



data <- plotCounts(dds,gene="SS1G_11992", intgroup=c("Hosts","TP"), 
                   returnData=TRUE,normalized = TRUE)

cond <- with(colData(dds), factor(paste(Hosts,TP)))
test_data <- assay(rld)[rownames(assay(rld))=="SS1G_05875",]
single_gene_rld <- as.data.frame(test_data)
single_gene_rld$group <- paste0(rld$Hosts,"-",rld$TP)
stripchart( single_gene_rld$test_data ~ single_gene_rld$group, method="jitter", 
            vertical=TRUE, las=2,pch=c(15,16,17),col=c("red","green","cyan"))
single_gene_rld <- as.data.frame(test_data)
single_gene_rld$Hosts <- paste0(rld$Hosts)
single_gene_rld$TP <- paste0(rld$TP)

test <- as.character(spep[spep$Gene_ID %in% "SS1G_12776",]$Annotation)
PICK48[rownames(PICK48)=="SS1G_05875",]
LTCK48[rownames(LTCK48)=="SS1G_05875",]
Trt_CK48[rownames(Trt_CK48)=="SS1G_05875",]


spep <- read.csv("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/SPEP.csv",header = TRUE,)
counts_normal <- as.data.frame(counts(dds,normalized=TRUE))
firstrow <- as.list(rownames(counts_normal))
GeneID_all <- unlist(firstrow)
counts_normal <- cbind(GeneID_all,counts_normal )
# 
# GeneID_all <- unlist(lapply(firstrow,function(x) {
#   ifelse(as.vector(unlist(strsplit(x,"\\|")))[2] %in% "Unknown",
#          as.vector(unlist(strsplit(x,"\\|")))[1],
#          as.vector(unlist(strsplit(x,"\\|")))[2])}))

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
draw_genename(T24PI,"")
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



#### Heatmap for top expressed genes ##########
topVarGenes <- head(order(-rowVars(assay(rld))),100)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
colors2 <- c("Blue","Yellow","Red")
sidecols <- c("grey","dodgerblue","lightgreen")[ rld$TP ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$Hosts,"-",rld$TP)
hc.rows1 <- hclust(dist(mat))
hc.cols1 <- hclust(dist(t(mat)))
heatmap.2(mat, trace="none", , ColSideColors=sidecols,col=colors,
           labRow=TRUE, key=FALSE,lhei=c(2,8), lwid=c(1,4),margins=c(6.5,4),
          Colv=as.dendrogram(hc.cols1))
#heatmap(mat,col=colors,cexRow = 0.5,ColSideColors = sidecols,)
#dev.off()
par(mfrow=c(1,1),cex = 1)
heatmap(mat[cutree(hc.rows1,k=6)==4,], Colv=as.dendrogram(hc.cols1), scale='none',
        col=colors,cexRow = 1,ColSideColors = sidecols,margins=c(7,7))
genenames <- rownames(mat[cutree(hc.rows1,k=6)==4,])
write.table(genenames,"clade-4",quote = FALSE,row.names = FALSE,col.names=FALSE)


##### Heat map for candidate effector genes ##########
spep_effector <- spep[spep$Effector_candidate %in% c("YES"),]
spep_effector_count <- data.frame()
for (prot in spep_effector$Gene_ID) {
  if (prot %in% counts_normal$GeneID_all) spep_effector_count <- rbind(spep_effector_count,counts_normal[GeneID_all %in% prot,])
}
spep_effector_count<- subset(spep_effector_count, select = -c(GeneID_all))
spep_effector_count <- as.matrix(spep_effector_count)
spep_effector_count <- spep_effector_count - rowMeans(spep_effector_count)
colnames(spep_effector_count) <- paste0(rld$Hosts,"-",rld$TP)
rownames(spep_effector_count)
top50 <- head(order(-rowVars(spep_effector_count)),50)
top_effector_count <- spep_effector_count[top50,]
heatmap.2(top_effector_count, trace="none", col=colors, ColSideColors=sidecols,
          lhei=c(.7,2), lwid=c(1.5,5), key=TRUE, keysize=0.02,revC=FALSE,margins=c(6.5,5.5),
          Colv=as.dendrogram(hc.cols1))
heatmap(top_effector_count,col=colors,scale = "none")
hc.rows2 <- hclust(dist(top_effector_count ))
hc.cols2 <- hclust(dist(t(top_effector_count )))
Effectors_sub <- rownames(top_effector_count[cutree(hc.rows2,k=2)==1,])

############ GO Enrichment ##############
rowsum.threshold <- 1 # user chosen
fdr.threshold <- 0.1 # user chosen
rs <- rowSums(counts(dds))
dds <- dds[ rs > rowsum.threshold ,]
dds <- DESeq(dds)
res <- results(dds, independentFiltering=FALSE) # use count threshold instead of IF
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ]
listMarts()
listDatasets(ensembl)


############ Donut chart to show gene categories in heatmap top100 DE genes ##############
# Input data
clade1_mole <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade1_pie_mole_funct.txt",
                          header=FALSE,sep="\t")
clade2_mole <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade2_pie_mole_funct.txt",
                          header=FALSE,sep="\t")
clade3_mole <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade3_pie_mole_funct.txt",
                          header=FALSE,sep="\t")
clade4_mole <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade4_pie_mole_funct.txt",
                          header=FALSE,sep="\t")
clade1_prot <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade1_pie_prot_class.txt",
                          header=FALSE,sep="\t")
clade2_prot <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade2_pie_prot_class.txt",
                          header=FALSE,sep="\t")
clade3_prot <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade3_pie_prot_class.txt",
                          header=FALSE,sep="\t")
clade4_prot <- read.table("/home/peascl/scl_de/htout_STAR_mergeGTF/R_output/data/clade4_pie_prot_class.txt",
                          header=FALSE,sep="\t")

par(mfrow=c(2,1),cex = 0.7)
draw_pie(clade1_mole,"Clade1")
draw_pie(clade2_mole,"Clade2")
draw_pie(clade3_mole,"Clade3")
draw_pie(clade4_mole,"Clade4")
par(mfrow=c(2,1),cex = 0.7)
draw_pie(clade1_prot,"Clade1")
draw_pie(clade2_prot,"Clade2")
draw_pie(clade3_prot,"Clade3")
draw_pie(clade4_prot,"Clade4")
