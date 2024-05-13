#######################################################
#Fig3 scenic
library(SingleCellExperiment)
library(SCENIC)

ace<-readRDS("EC_sce.rds")

exprMat <- counts(ace)

cellInfo <- colData(ace)

cellInfo<-data.frame(row.names=rownames(cellInfo),CellType=cellInfo$SubCluster)

saveRDS(cellInfo, file="int/cellInfo.Rds")

exprMat<-as.matrix(exprMat)

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(org="hgnc", dbDir="/home/chuxj/scenic", nCores=10)

scenicOptions@inputDatasetInfo$cellInfo<-"int/cellInfo.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions)

exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#######################################################
##cell subset co-occurence
setwd("/data1/chuxj/project/allCells/")

library(Seurat)
library(ggplot2)


dat<-readRDS("AllAnnotated.RDS")

dat<-subset(dat,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib","Glial","ILC","MAST",
                                       "Mono/Macro","Plasma","Proliferating Myeloids","NK","Proliferating T"))

dat1<-subset(dat,Dataset %in% c("C.Smillie","G.Li","GSE132257","GSE132465",
                                "GSE144735","GSE150115","GSE188711","GSE201349",
                                "J.Qi","J.Qian","Pelka.K","R.Elmentaite"))

tmp<-table(dat1$Sample,dat1$SubCluster)

Pa<-which(apply(tmp,1,sum)>200)
tmp<-tmp[Pa,]
tmp<-tmp/apply(tmp,1,sum)

meta<-read.table("../tmp/sample_anno.txt",header=T,row.names = 1,sep="\t")

library(Hmisc)

tmp<-tmp[,c(1:21,24:29,31:56)]

write.table(tmp,file="Sample_prop.txt",quote=F,sep="\t")

tmp<-read.table("Sample_prop.txt",header=T,row.names = 1,sep="\t")

cor<-rcorr(as.matrix(tmp),type = "spearman")

library(pheatmap)

col=c(rev(brewer.pal(9,"Blues")),"White",brewer.pal(9,"Reds"))

p0=pheatmap(cor$r,breaks = seq(-1,1,0.1),color = col,clustering_method = "ward.D")

#####
meta<-meta[rownames(tmp),]
Imm1<-tmp[which(meta$Class=="T"),]
Imm2<-tmp[which(meta$Class=="N"),]

####
Imm_cor1<-rcorr(as.matrix(Imm1),type = "spearman")
Imm_cor2<-rcorr(as.matrix(Imm2),type = "spearman")

col=c(rev(brewer.pal(9,"Blues")),"White",brewer.pal(9,"Reds"))

p1=pheatmap(Imm_cor1$r,breaks = seq(-1,1,0.1),color = col,clustering_method = "ward.D")
p2=pheatmap(Imm_cor2$r,breaks = seq(-1,1,0.1),color = col,cluster_rows = p1$tree_row,cluster_cols = p1$tree_col)

