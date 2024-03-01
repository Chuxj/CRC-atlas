#############################
#Fig6 de analysis

library(Seurat)
Sobj<-readRDS("AllAnnotated.RDS")
Sobj<-subset(Sobj,Class %in% c("N","T"))
Idents(Sobj)<-"Class"
Sobj<-subset(Sobj,ParentalCluster %in% c("B","CD4","CD8","DC","EC","Fib","Glial","ILC",
                                         "MAST","Mono/Macro","NK","Plasma"))

seed=1
set.seed(seed)
  
  for (i in c("B","CD4","CD8","DC","EC","Fib","Glial","ILC",
              "MAST","Mono/Macro","NK","Plasma")){
    
    message(paste0("working with ",i))
    
    tmp1 <- colnames(Sobj)[which(Sobj$ParentalCluster ==i & Sobj$Class=="N")]
    tmp2 <- colnames(Sobj)[which(Sobj$ParentalCluster ==i & Sobj$Class=="T")]
    
    if(length(tmp1)>1000){cells1<-sample(tmp1,1000)}else{cells1<-tmp1}
    if(length(tmp2)>1000){cells2<-sample(tmp2,1000)}else{cells2<-tmp2}
    
    markers <- FindMarkers(subset(Sobj,cells=c(cells1,cells2)), ident.1 = "T", ident.2 = "N", min.pct = 0.1,logfc.threshold=0)  
    
    if (i == "Mono/Macro"){
      write.table(markers, file=paste0("/data1/chuxj/project/gwas/seed",seed,"/DE_TvN_MoMac.txt"),quote=F,sep="\t")
    }else{
      write.table(markers, file=paste0("/data1/chuxj/project/gwas/seed",seed,"/DE_TvN_",i,".txt"),quote=F,sep="\t")}
  }

#############################
#qqplot
dat<-read.table("Stat.txt",header=T,sep = "\t")

EXPPVAL<-function(my.pvalues){(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)}

p_group<-split(dat,dat$CellType)

y<-lapply(p_group,function(w){p<-w$Pval;EXPPVAL(p)})

p<-c(y$B,y$CD4,y$CD8,y$DC,y$EC,y$Fib,y$Glial,y$ILC,y$MAST,y$MoMac,y$NK,y$Plasma)
Subpopulation<-c(rep('B',length(y$B)),rep('CD4',length(y$CD4)),rep("CD8",length(y$CD8)),
                 rep("DC",length(y$DC)),rep("EC",length(y$EC)),rep("Fib",length(y$Fib)),
                 rep("Glial",length(y$Glial)),rep("ILC",length(y$ILC)),
                 rep("MAST",length(y$MAST)),rep("MoMac",length(y$MoMac)),
                 rep("NK",length(y$NK)),rep("Plasma",length(y$Plasma))
)
exp<-data.frame(p,Subpopulation)
exp<-transform(exp,logp=(-log10(p)))

or<-order(dat$CellType)
dat<-dat[or,]

dat<-transform(dat,logp=(-log10(Pval)))

tmp<-data.frame('obs'=dat$logp,'exp'=exp$logp,CellType=dat$CellType)

tmp$CellType<-factor(tmp$CellType,levels=c("B","Plasma","CD4","CD8","NK","ILC","DC","MoMac","MAST","EC","Fib","Glial"))

mycol_CellType=c("#e41a1c","#e7298a","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#e5c494","#66c2a5","#f781bf","#999999")

ggplot(data=tmp,aes(x=exp,y=obs,color=CellType))+geom_point()+geom_abline(slope = 1,intercept = 0,lty=2,color="grey20")+
  theme_classic()+scale_color_manual(values = mycol_CellType)

#############################
#lambda

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

