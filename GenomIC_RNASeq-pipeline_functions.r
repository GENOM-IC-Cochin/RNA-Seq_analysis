normalization.MatrixFile <- function(dds,ensembl.database,ensembl.attributes,projectName){
  
  countNorm<-counts(dds, normalized=TRUE)
  ensemblIDs <- row.names(countNorm)
  
  IDsWithNamesDesc <- getBM(attributes = ensembl.attributes, 
                            filters = 'ensembl_gene_id', 
                            values = ensemblIDs, 
                            mart = ensembl.database,
                            useCache = FALSE)
  IDsWithNamesDesc %<>% mutate(across(everything(), na_if, ""))
  rownames(IDsWithNamesDesc) <- make.names(IDsWithNamesDesc$ensembl_gene_id, unique=TRUE)
  dataMerged <- merge(as.data.frame(countNorm),IDsWithNamesDesc,by="row.names",all.x=TRUE)
  nameNorm <- paste(projectName,"_deseq2_NormalizedMatrix.tsv",sep="")
  write.table(dataMerged,file=nameNorm,sep='\t',row.names=F)
  
  return(IDsWithNamesDesc)
  
}

normalization.boxPlot <- function (countNorm,projectName){
  pseudoNormCount <- as.data.frame(log2(countNorm+1),row.names=row.names(countNorm))
  pseudoNormCount$Ids <- row.names(pseudoNormCount)
  datNorm <- melt(pseudoNormCount, id.vars = "Ids", variable.name = "Samples", value.name = "count")
  colnames(datNorm)<-c("Ids", "Samples", "count")
  datNorm <- merge(datNorm, configuration, by.x = "Samples", by.y = "Name")
  
  countRaw<-counts(dds, normalized=FALSE)
  pseudoRawCount <- as.data.frame(log2(countRaw+1),row.names=row.names(countRaw))
  pseudoRawCount$Ids <- row.names(pseudoRawCount)
  datRaw <- melt(pseudoRawCount, id.vars = "Ids", variable.name = "Samples", value.name = "count")
  colnames(datRaw)<-c("Ids", "Samples", "count")
  datRaw <- merge(datRaw, configuration, by.x = "Samples", by.y = "Name")
  
  gDens <- ggplot(datNorm, aes(x = count, colour = Samples, fill = Samples)) + geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) + theme(legend.position = "top") + xlab(expression(log[2](count + 1)))
  ggsave(filename=paste(projectName,"_DensityNormCount.jpg",sep=""), plot=gDens)
  
  gNorm.name <- paste(projectName,"_boxplotNormCount.jpg",sep="")
  gNorm <- ggplot(data = datNorm, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + ggtitle("Normalized log2(counts) per sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename=gNorm.name, plot=gNorm)
  
  gRaw.name <- paste(projectName,"_boxplotRawCount.jpg",sep="")
  gRaw <- ggplot(data = datRaw, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + ggtitle("Raw log2(counts) per sample") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename=gRaw.name, plot=gRaw)
  
  return(c(gRaw.name,gNorm.name))
}

unsupervised.PCA <- function (rld,dds,projectName){
  ntop <- 500
  
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(rld)[select, ] )
  pc <- prcomp(mat)
  eig <- (pc$sdev)^2
  variance <- eig*100/sum(eig)
  
  PCAdata<-as.data.frame(pc$x[,1:3])
  PCAdata$condition <- dds$Condition
  
  nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes.jpg",sep="")
  pca1 <- ggplot(PCAdata,aes(x=PC1,y=PC2,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC2: ",round(variance[2],1),"% variance"))
  ggsave(filename=nameACP1, plot=pca1)
  
  nameACP2 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC2vsPC3_top500varGenes.jpg",sep="")
  pca2 <- ggplot(PCAdata,aes(x=PC2,y=PC3,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC2: ",round(variance[2],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP2, plot=pca2)
  
  nameACP3 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC3_top500varGenes.jpg",sep="")
  pca3 <- ggplot(PCAdata,aes(x=PC1,y=PC3,label=rownames(PCAdata))) + geom_point(aes(shape=condition, color=condition), size = 5) + geom_point() + geom_label_repel() + xlab(paste0("PC1: ",round(variance[1],1),"% variance")) + ylab(paste0("PC3: ",round(variance[3],1),"% variance"))
  ggsave(filename=nameACP3, plot=pca3)
  
  # nameACP_contrib <- paste(projectName,"_deseq2_Unsupervised_PCA_contribution.png",sep="")
  # contrib <- fviz_eig(pc)
  # ggsave(filename=nameACP_contrib, plot=contrib)
  
  return(c(nameACP1,nameACP2,nameACP3))
}

unsupervised.HC <- function (rld,configuration,projectName){
  distsRL <- dist(t(assay(rld)))
  matDist <- as.matrix(distsRL)
  hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
  Conditions <- data.frame(configuration$Condition,row.names=configuration$Name)
  nameClustering <- paste(projectName,"deseq2_Unsupervised_clustering_euclidean-complete.png",sep="_")
  png(filename=nameClustering,width=7 ,height=7, units="in",res = 600 ) 
  pheatmap(matDist, 
           col=hmcol, 
           annotation_col = Conditions,
           show_rownames=F,
           show_colnames=T)
  dev.off()
  
  nameClustering <- paste(projectName,"deseq2_Unsupervised_clustering_euclidean-complete_no-sampleNames.png",sep="_")
  png(filename=nameClustering,width=7 ,height=7, units="in",res = 600 ) 
  pheatmap(matDist, 
           col=hmcol, 
           annotation_col =  Conditions,
           show_rownames=F,
           show_colnames=F)
  dev.off()
}

supervised.res <- function (contrasteList,dds,projectName,IDsWithNamesDesc,count.normalized,all_results,z,list_volcanoPlot){
  ctrst <- unname(unlist(contrasteList[z,]))
  res <- results(dds, contrast=ctrst,alpha=0.05)
  name <- paste(projectName,"deseq2_results_contrast",ctrst[2],'vs',ctrst[3],sep="_")
  res <- merge(as.data.frame(res),IDsWithNamesDesc,by="row.names",all.x=TRUE)
  row.names(res)<-res$Row.names
  res<-res[,2:10]
  fullData <- merge(as.data.frame(res),count.normalized,by="row.names",all.x=TRUE)
  row.names(fullData)<-fullData$Row.names
  fullData<-fullData[order(fullData$padj), ]
  all_results[[z]] <- fullData
  #names(all_results)[z] <- paste(ctrst[2], "vs", ctrst[3], sep = "_")
  total_tests <- nrow(fullData)
  fullData$pdaj.bonferroni <- p.adjust(fullData$pvalue, method = "bonferroni", n = total_tests)
  write.table(as.data.frame(fullData),file=paste(name,".tsv",sep=""),sep='\t',row.names=F)
  volcano_plot <- make_nice_volcanoPlot(fullData,name,ctrst)
  make_nice_diffPlot(dds,fullData,configuration,name,ctrst) 
  list_volcanoPlot <- c(list_volcanoPlot,volcano_plot)
  return(all_results)
}

make_nice_volcanoPlot <- function(fullData,name,ctrst) {
  
  
  newdata<-fullData[,c(1,3,7,9)]
  
  volcanoTitle <- paste(ctrst[2],'vs',ctrst[3],sep=" ")
  
  colnames(newdata)<-c("Id","log2FoldChange","padj","GeneSymbol")
  
  # add a column of NAs
  newdata$DEGs <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  newdata$DEGs[newdata$log2FoldChange > log2(1.5)] <- "UP"
  newdata$DEGs[newdata$log2FoldChange > log2(1.5) & newdata$padj < 0.05] <- "UP & DE"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  newdata$DEGs[newdata$log2FoldChange < -log2(1.5)] <- "DOWN"
  newdata$DEGs[newdata$log2FoldChange < -log2(1.5) & newdata$padj < 0.05] <- "DOWN & DE"
  
  nameVolcanoPlot<-paste(name,"_volcanoPlot.png",sep="")
  
  volcanoNoLabel<-ggplot(newdata) +
    aes(x = log2FoldChange, y = -log10(padj), colour = DEGs) +
    geom_point(shape = "circle", size = 1.5) +
    xlab("log2(Fold-Change") +
    ylab("-log10(padj)") +
    labs(title = volcanoTitle) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="grey50",linetype=3) + 
    geom_hline(yintercept=-log10(0.05), col="grey50",linetype=3) +
    scale_color_manual(values=c("#CCCCFF","blue3", "snow3", "#FFCCCC", "firebrick2"))
  
  ggsave(filename=nameVolcanoPlot, plot=volcanoNoLabel)
  
  top10DEG<-newdata[1:10,]
  
  nameVolcanoPlot<-paste(name,"_volcanoPlot_Top10.png",sep="")
  
  volcanoTop10<-ggplot(newdata) +
    aes(x = log2FoldChange, y = -log10(padj), colour = DEGs) +
    geom_point(shape = "circle", size = 1.5) +
    xlab("log2(Fold-Change") +
    ylab("-log10(padj)") +
    labs(title = volcanoTitle ,subtitle = "Top 10 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="grey50",linetype=3) + 
    geom_hline(yintercept=-log10(0.05), col="grey50",linetype=3) +
    scale_color_manual(values=c("#CCCCFF","blue3", "snow3", "#FFCCCC", "firebrick2"))+
    geom_label_repel(data=top10DEG, aes(label=GeneSymbol), color = 'black',show.legend=FALSE)
  ggsave(filename=nameVolcanoPlot, plot=volcanoTop10)
  
  
}

make_nice_diffPlot <- function(dds,fullData,configuration,name,ctrst) {
  
  diffTitle <- paste(ctrst[2],'vs',ctrst[3],sep=" ")
  
  countNorm<-counts(dds, normalized=TRUE)
  pseudoNormCount <- as.data.frame(log2(countNorm+1),row.names=row.names(countNorm))
  pseudoNormCount$Ids <- row.names(pseudoNormCount)
  
  ListTop_10 <- fullData[1:10,1]
  
  Top10DEGs <- subset(countNorm, row.names(countNorm) %in% ListTop_10)
  Top10DEGs <- melt(Top10DEGs)
  colnames(Top10DEGs) <- c("EnsemblID","Sample","NormalizedCount")
  
  top_10 <- merge(Top10DEGs, configuration, by.x = "Sample",by.y="Name")
  
  top_10 <- subset(top_10, Condition %in% ctrst)
  
  diffplot <- ggplot(top_10) +
    geom_point(aes(x = EnsemblID, y = NormalizedCount, color = Condition, shape=Condition)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("Normalized Counts") +
    labs(title = diffTitle ,subtitle = "Top 10 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  nameDiffPlot<-paste(name,"_diffPlot_Top10.png",sep="")
  ggsave(filename=nameDiffPlot, plot=diffplot)
}
