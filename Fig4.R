###################################################
#Fig4 classification

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
#Fig4 cellchat
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

###########
##EMT and TGFB signaling score

setwd("/data1/chuxj/project/otherClassification")

count<-read.table("psudobulk.txt",header=T,row.names = 1,sep="\t")

anno_pat<-read.table("/data1/chuxj/project/cellModule/classification/Patient_Group6.txt",header=T,row.names = 1,sep="\t",check.names = F)

ol<-intersect(names(count),rownames(anno_pat))

count<-count[,ol]
anno_pat<-anno_pat[ol,]
anno_pat<-as.data.frame(anno_pat)
names(anno_pat)<-"Group"
rownames(anno_pat)<-ol

library(Seurat)
Sobj<-CreateSeuratObject(count,meta.data = anno_pat)
Sobj<-NormalizeData(Sobj)

GO_EMT<-read.table("../classification/EMT_TGFB/GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt")

GO_TGFB<-read.table("../classification/EMT_TGFB/GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION.txt")

######
glist<-list(GO_EMT=GO_EMT$V1,GO_TGFB=GO_TGFB$V1)

Sobj<-AddModuleScore(Sobj,glist)

tmp<-FetchData(Sobj,vars = c("Cluster1","Cluster2","Group"))

names(tmp)<-c("GO_EMT","GO_TGFB","Group")

tmp$Group<-factor(tmp$Group,levels = c(1,2,3,4,5,6))

comparisons=list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6))

p1=ggplot(tmp,aes(x=Group,y=GO_EMT,color=Group))+geom_boxplot()+theme_classic()+geom_jitter(size=0.5)+scale_color_manual(values = c("#8f5362","#d37b6d","#9c8207","#f6df56","#79b4a0","#706d94"))+
  stat_compare_means(comparisons = comparisons,method = "t.test")

p2=ggplot(tmp,aes(x=Group,y=GO_TGFB,color=Group))+geom_boxplot()+theme_classic()+geom_jitter(size=0.5)+scale_color_manual(values = c("#8f5362","#d37b6d","#9c8207","#f6df56","#79b4a0","#706d94"))+
  stat_compare_means(comparisons = comparisons,method = "t.test")
