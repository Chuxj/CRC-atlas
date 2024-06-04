##################################
#Fig5 MeanZscore

Sobject<-readRDS("TumorSample_add6Group.RDS")

Sobject<-subset(Sobject,Group %in% c(1,2,3,4,5,6))

Idents(Sobject)<-"Group"

library(pheatmap)
library(RColorBrewer)

col<-c(brewer.pal(9,"RdBu"))
col=colorRampPalette(col)(100)

#####CD8
gene=c("BTLA","TNFRSF18","CTLA4","TIGIT","PDCD1","HAVCR2","LAG3")

mat<-FetchData(subset(Sobject,ParentalCluster=="CD8"),vars = c(gene,"Group"))
test<-apply(mat[,-ncol(mat)],2,scale)
rownames(test)<-rownames(mat)
test<-as.data.frame(test)
x<-split(test,mat$Group)
x2<-lapply(x, function(w){apply(w,2,mean)})
x<-do.call(rbind,x2)

pheatmap(t(x),cluster_rows = F,cluster_cols = F,color=rev(col),breaks=seq(-1.8,1.8,0.036))

#####malignant
gene=c("CD274","PDCD1LG2","CD47")

mat<-FetchData(subset(Sobject,SubCluster=="Malignant cells"),vars = c(gene,"Group"))
test<-apply(mat[,-ncol(mat)],2,scale)
rownames(test)<-rownames(mat)
test<-as.data.frame(test)
x<-split(test,mat$Group)
x2<-lapply(x, function(w){apply(w,2,mean)})
x<-do.call(rbind,x2)

pheatmap(t(x),cluster_rows = F,cluster_cols = F,color=rev(col),breaks=seq(-1.8,1.8,0.036))

#####Macroph
gene=c("SIRPA","FCGR1A","FCGR2A","FCGR2B","FCGR3A","FCGR3B")

mat<-FetchData(subset(Sobject,SubCluster %in% c("Macro-C1QC","Macro-ISG15","Macro-LYVE1",
                                                "Macro-SPP1")),vars = c(gene,"Group"))
test<-apply(mat[,-ncol(mat)],2,scale)
rownames(test)<-rownames(mat)
test<-as.data.frame(test)
x<-split(test,mat$Group)
x2<-lapply(x, function(w){apply(w,2,mean)})
x<-do.call(rbind,x2)

pheatmap(t(x),cluster_rows = F,cluster_cols = F,color=rev(col),breaks=seq(-1.8,1.8,0.036))

##################################
#Fig5 CD8-CXCL13, CD4-CXCL13

library(Seurat)

Sobj<-readRDS("AllAnnotated.RDS")

c1<-subset(Sobj,ParentalCluster=="CD8")
c2<-subset(Sobj,ParentalCluster=="CD4")

Idents(c1)<-"SubCluster"
Idents(c2)<-"SubCluster"

c1.markers <- FindMarkers(c1, ident.1 = "CD8-CXCL13", ident.2 = "CD8-IL7R", min.pct = 0.25,logfc.threshold = 0)
c2.markers <- FindMarkers(c2, ident.1 = "CD4-CXCL13", ident.2 = "CD4-CCR7", min.pct = 0.25,logfc.threshold = 0)

ol<-intersect(rownames(c1.markers),rownames(c2.markers))

plot<-data.frame(c1=c1.markers[ol,"avg_log2FC"],c2=c2.markers[ol,"avg_log2FC"],
                 c1p=c1.markers[ol,"p_val_adj"],c2p=c2.markers[ol,"p_val_adj"])
rownames(plot)<-ol

library(ggplot2)
library(ggrepel)

plot<-transform(plot,gene=rownames(plot))

plot$add<-"No"
plot$add[which(plot$gene=="TIGIT")]<-"Yes"
plot$add[which(plot$gene=="PDCD1")]<-"Yes"
plot$add[which(plot$gene=="CXCL13")]<-"Yes"
plot$add[which(plot$gene=="TNFRSF18")]<-"Yes"
plot$add[which(plot$gene=="CTLA4")]<-"Yes"

options(ggrepel.max.overlaps = Inf)

ggplot(plot,aes(x=c1,y=c2,color=add))+geom_point(alpha=0.6)+theme_classic()+ylim(-2.5,4)+xlim(-2.5,4)+
  geom_hline(yintercept = 0,lty="dashed",color="grey56")+
  geom_vline(xintercept = 0,lty="dashed",color="grey56")+
  geom_smooth(method='lm',color="#bd706e")+
  geom_text_repel(aes(label = ifelse(add  == "Yes", as.character(gene), '')))+
  scale_colour_manual(name = "add", 
                      values = c("Yes" = "#b2182b","No"="#2166ac"))+
  theme(legend.position = "None")+
  ylab("Average LFC in CD4-CXCL13") +
  xlab("Average LFC in CD8-CXCL13")

cor.test(plot$c1,plot$c2)
