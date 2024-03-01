###################################################
#Fig3 classification

library(pheatmap)
library(RColorBrewer)

col=c("white",brewer.pal(9,"Reds"))
col<-c(brewer.pal(9,"RdBu"))

col_test<-c(col[1:2],"white",col[8:9])

pat_meta<-read.table("/data1/chuxj/project/cellAnno/pat_fullAnno.txt",sep="\t",header=T,row.names = 1)


pat_meta<-pat_meta[which(pat_meta$Dataset%in% c("GSE132257","GSE132465","GSE144735","GSE188711"
                                                ,"GSE201349","J.Qi","J.Qian","Pelka.K")),]

mat<-read.table("../cellModule/classification/allcells/Patient_proportion_matrix_allcells.txt",header=T,row.names = 1,sep="\t")


ol<-intersect(names(mat),rownames(pat_meta))

mat<-mat[,ol]

anno<-data.frame(CellType=rep(NA,nrow(mat)),MajorCellType=rep(NA,nrow(mat)),row.names =rownames(mat) )
anno$CellType[1:5]<-"CD4 T"
anno$CellType[6:11]<-"CD8 T"
anno$CellType[13:14]<-"NK"
anno$CellType[16:21]<-"B"
anno$CellType[22:33]<-"Myeloid cells"
anno$MajorCellType[1:33]<-"CD45pos"

anno$CellType[34:45]<-"EC"
anno$CellType[46:55]<-"Fib"
anno$MajorCellType[34:56]<-"CD45neg"

mycol_MajorCellType=c("#3b6799","#d5ba91")
mycol_CellType=brewer.pal(n=8,"Dark2")

names(mycol_MajorCellType)<-as.character(unique(anno$MajorCellType))
names(mycol_CellType)<-as.character(unique(anno$CellType))

annotation_colors<-list(MajorCellType=mycol_MajorCellType,CellType=mycol_CellType)

mat<-mat[c(1:11,13:17,19:32,34:37,40:44,46:50,53:56),]
anno<-anno[c(1:11,13:17,19:32,34:37,40:44,46:50,53:56),]

col=colorRampPalette(rev(col_test))(50)

p=pheatmap::pheatmap(as.matrix(mat),clustering_method = "ward.D",scale = "row",color = col
                     ,show_colnames = F,annotation_row = anno)

###################################################
#Fig3 cellchat
library(dplyr)
library(Seurat)
library(patchwork)
library(CellChat)

Sobject<-readRDS("TumorSample_add6Group.RDS")

Sobject<-subset(Sobject,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib"
                                               ,"Malignant cells","Mono/Macro","NK","Plasma"))

for (i in c(1,2,3,4,5,6)){
  
  message(paste0("Processing G",i))
  
  tmp<-subset(Sobject,Group == i)
  
  data.input <- GetAssayData(tmp, assay = "RNA", slot = "data") # normalized data matrix
  labels <- tmp$ParentalCluster
  meta <- data.frame(CellType = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
  
  CellChatDB <-CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
##note! use a low threshold here, and filter results afterwards
  
  df.net <- subsetCommunication(cellchat)
  
}
