#####################################################
#Fig1 UMAP

library(Seurat)
library(RColorBrewer)
library(ggplot2)

Sobj<-readRDS("AllAnnotated.RDS")

col2=c("#81C0C0","#FF6600","#02C874","#B766AD","#004897","#FF99CC","#FFD306")
DimPlot(Sobj, group.by = "Class",shuffle = T)+scale_color_manual(values=col2)

Idents(Sobj)<-"SubCluster"

col<-c(brewer.pal(12,"Set3"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"),brewer.pal(9,"Set1"),brewer.pal(8,"Dark2"),
       brewer.pal(8,"Pastel2"),brewer.pal(8,"Set2"),brewer.pal(8,"Paired"))

DimPlot(Sobj,shuffle = T,group.by = "SubCluster")+scale_color_manual(values=col)

#####################################################
#Fig1 boxplots

library(Seurat)
library(ggplot2)


dat<-readRDS("AllAnnotated.RDS")

dat<-subset(dat,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib","Glial","ILC","MAST",
                                       "Mono/Macro","Plasma","Proliferating Myeloids","NK","Proliferating T"))

dat1<-subset(dat,Dataset %in% c("C.Smillie","G.Li","GSE132257","GSE132465",
                                "GSE144735","GSE150115","GSE188711","GSE201349",
                                "J.Qi","J.Qian","Pelka.K","R.Elmentaite"))

tmp<-table(dat1$Sample,dat1$ParentalCluster)

Pa<-which(apply(tmp,1,sum)>200)
tmp<-tmp[Pa,]
tmp<-tmp/apply(tmp,1,sum)

Class<-data.frame(Sample=dat1$Sample,Class=dat1$Class)
flag<-duplicated(Class$Sample)
Class<-Class[which(flag=="FALSE"),]
rownames(Class)<-Class$Sample

plot<-transform(tmp,Class=Class[rownames(tmp),"Class"])

plot$Class<-factor(plot$Class,levels = c("Healthy","Uninflamed","Inflamed","N","Polyp","T"))
col_class<-c("#81C0C0","#FF99CC","#FF6600","#02C874","#B766AD","#004897")

ggplot(data=plot,aes(x=Class,y=Freq,fill=Class))+geom_boxplot()+geom_point(size=0.8)+scale_fill_manual(values=col_class)+theme_classic()+facet_grid(.~Var2,scales = "free")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

###
dat<-readRDS("AllAnnotated.RDS")
table(dat$ParentalCluster)
dat<-subset(dat,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib","Glial","ILC","MAST",
                                       "Mono/Macro","Plasma","Proliferating Myeloids","NK","Proliferating T"))

dat1<-subset(dat,Dataset %in% c("C.Smillie","G.Li","GSE132257","GSE132465",
                                "GSE144735","GSE150115","GSE188711","GSE201349",
                                "J.Qi","J.Qian","Pelka.K","R.Elmentaite"))

tmp<-table(dat1$Sample,dat1$SubCluster)

Pa<-which(apply(tmp,1,sum)>200)
tmp<-tmp[Pa,]
tmp<-tmp/apply(tmp,1,sum)

Class<-data.frame(Sample=dat1$Sample,Class=dat1$Class)
flag<-duplicated(Class$Sample)
Class<-Class[which(flag=="FALSE"),]
rownames(Class)<-Class$Sample

plot<-transform(tmp,Class=Class[rownames(tmp),"Class"])

plot$Class<-factor(plot$Class,levels = c("Healthy","Uninflamed","Inflamed","N","Polyp","T"))
col_class<-c("#81C0C0","#FF99CC","#FF6600","#02C874","#B766AD","#004897")

plotneg<-plot[which(plot$Var2 %in% c("Fib-CCL19","Fib-FABP5","Fib-IGF1","Fib-MFAP5","Fib-POSTN","Fib-WNT5A",
                                     "Pericyte","Smooth Muscle","Glial cells",
                                     "Artery-FBLN2","Artery-GJA5","Capillary-BTNL9","Capillary-CA4","HEV-CXCL10","HEV-SELE",
                                     "Lymmphatic ECs-LYVE1","Lymmphatic ECs-LowFlow","Vein-ACKR1")),]

plotpos<-plot[which(plot$Var2 %in% c("B-IgD","B-LRMP","B-MS4A1","Plasma-IgA","Plasma-IgG",
                                     "CD4-ANXA1","CD4-CCR7","CD4-CXCL13","CD4-IL17A","Treg-FOXP3",
                                     "CD8-CXCL13","CD8-gdT","CD8-GZMK","CD8-IL7R","CD8-ISG15","CD8-MAIT",
                                     "NK-GZMH","NK-XCL1","ILCs",
                                     "cDC-CD1C","cDC-CLEC9A","cDC-LAMP3","pDC","Mono-CD16","Mono-FCN1",
                                     "Macro-C1QC","Macro-ISG15","Macro-LYVE1","Macro-SPP1","MAST cells")),]

plotneg$Var2<-factor(plotneg$Var2,levels=c("Fib-CCL19","Fib-FABP5","Fib-IGF1","Fib-MFAP5","Fib-POSTN","Fib-WNT5A",
                                           "Pericyte","Smooth Muscle","Glial cells",
                                           "Artery-FBLN2","Artery-GJA5","Capillary-BTNL9","Capillary-CA4","HEV-CXCL10","HEV-SELE",
                                           "Lymmphatic ECs-LYVE1","Lymmphatic ECs-LowFlow","Vein-ACKR1"))

plotpos$Var2<-factor(plotpos$Var2,levels=c("B-IgD","B-LRMP","B-MS4A1","Plasma-IgA","Plasma-IgG",
                                           "CD4-ANXA1","CD4-CCR7","CD4-CXCL13","CD4-IL17A","Treg-FOXP3",
                                           "CD8-CXCL13","CD8-gdT","CD8-GZMK","CD8-IL7R","CD8-ISG15","CD8-MAIT",
                                           "NK-GZMH","NK-XCL1","ILCs",
                                           "cDC-CD1C","cDC-CLEC9A","cDC-LAMP3","pDC","Mono-CD16","Mono-FCN1",
                                           "Macro-C1QC","Macro-ISG15","Macro-LYVE1","Macro-SPP1","MAST cells"))

ggplot(data=plotpos,aes(x=Class,y=Freq,fill=Class))+geom_boxplot()+geom_point(size=0.8)+scale_fill_manual(values=col_class)+theme_classic()+facet_grid(.~Var2,scales = "free")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

ggplot(data=plotneg,aes(x=Class,y=Freq,fill=Class))+geom_boxplot()+geom_point(size=0.8)+scale_fill_manual(values=col_class)+theme_classic()+facet_grid(.~Var2,scales = "free")+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())

#####################################################
#pca

library(Seurat)
library(ggplot2)

setwd("/data1/chuxj/project/allCells")

dat<-readRDS("AllAnnotated.RDS")

dat<-subset(dat,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib","Glial","ILC","MAST",
                                       "Mono/Macro","Plasma","Proliferating Myeloids","NK","Proliferating T"))

dat1<-subset(dat,Dataset %in% c("C.Smillie","G.Li","GSE132257","GSE132465",
                                "GSE144735","GSE150115","GSE188711","GSE201349",
                                "J.Qi","J.Qian","Pelka.K","R.Elmentaite"))

dat1$Class.Patient=paste(dat1$Class,dat1$Patient,sep = "-")
tmp<-table(dat1$Class.Patient,dat1$SubCluster)
tmp<-tmp/apply(tmp,1,sum)

mat<-prcomp(tmp)

pcVar<-summary(mat)
varPC1<-pcVar$importance[2,1]
varPC2<-pcVar$importance[2,2]

mat<-mat$x[,1:2]

mat<-as.data.frame(mat)

anno<-read.table("/data1/chuxj/project/cellAnno/pat_class_anno.txt",header=T,row.names = 1,sep="\t")

anno$Site<-factor(anno$Site,levels=c("left","right"))

ol<-intersect(rownames(anno),rownames(mat))

mat$Site<-anno[ol,"Site"]
mat$Dataset<-anno[ol,"Dataset"]
mat$Class<-anno[ol,"Class"]

col_class<-c("#00ffff","#FF99CC","#FF6600","#02C874","#B766AD","#004897")
col_Site<-c("#00bbff","#ffb7dd")
col_Dataset=c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"))

mat$Class<-factor(mat$Class,levels=c("Healthy","Uninflamed","Inflamed","N","Polyp","T"))

p1=ggplot(data=mat,aes(x=PC1,y=PC2,color=Class))+geom_point()+scale_color_manual(values = col_class)+theme_classic()+
  xlab(paste("PC1 explained variance=",varPC1*100,"%",sep=""))+ylab(paste("PC2 explained variance=",varPC2*100,"%",sep=""))

p2=ggplot(data=mat,aes(x=PC1,y=PC2,color=Site))+geom_point()+scale_color_manual(values = col_Site)+theme_classic()+
  xlab(paste("PC1 explained variance=",varPC1*100,"%",sep=""))+ylab(paste("PC2 explained variance=",varPC2*100,"%",sep=""))

p3=ggplot(data=mat,aes(x=PC1,y=PC2,color=Dataset))+geom_point()+scale_color_manual(values = col_Dataset)+theme_classic()+
  xlab(paste("PC1 explained variance=",varPC1*100,"%",sep=""))+ylab(paste("PC2 explained variance=",varPC2*100,"%",sep=""))

#distance
test<-split(mat[,1:2],mat$Class)
test<-lapply(test,function(w){apply(w,2,mean)})
test<-do.call(rbind,test)

d<-dist(test)

d<-as.matrix(d)

library(corrplot)

corrplot(d,type="lower",is.corr = FALSE,col.lim=c(0,0.5),addCoef.col="grey55",number.font = 0.05,method = "color")
