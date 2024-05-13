#####################################################
# Fig2 bh dist
library(distdimscr)
library(tidyverse)

Sobj<-readRDS("AllAnnotated.RDS")

#remove dataset with only CD45+
Sobj<-subset(Sobj,Dataset!="GSE125527")

### T vs. Inf
cellType<-c("B","Plasma","CD4","CD8","Fib","Mono/Macro","NK","EC")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){
  
  for (j in 1:10){
    
    set.seed(j)
    
    n=length(which(Sobj$Class=="Inflamed" & Sobj$ParentalCluster==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Class=="T" & Sobj$ParentalCluster==CT)],n)
    
    cells.Inf <- colnames(Sobj)[Sobj$ParentalCluster==CT & Sobj$Class=="Inflamed"]
    
    tmp<-Sobj@reductions$harmony@cell.embeddings
    
    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]
    
    for (i in 1:10) {
      
      d<-(j-1)*10+i
      
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)
      
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)
      
    }
  }
}

### T vs. Polyp

cellType<-c("B","CD4","CD8","Fib","Mono/Macro")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){
  
  for (j in 1:10){
    
    set.seed(j)
    
    n=length(which(Sobj$Class=="Polyp" & Sobj$ParentalCluster==CT))
    cells.T<-sample(colnames(Sobj)[which(Sobj$Class=="T" & Sobj$ParentalCluster==CT)],n)
    
    cells.Inf <- colnames(Sobj)[Sobj$ParentalCluster==CT & Sobj$Class=="Polyp"]
    
    tmp<-Sobj@reductions$harmony@cell.embeddings
    
    cells.Inf.pca <- tmp[cells.Inf,]
    cells.T.pca <- tmp[cells.T,]
    
    for (i in 1:10) {
      
      d<-(j-1)*10+i
      
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)
      
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)
      
    }
  }
}

#### Para vs. T
cellType<-c("B","CD4","CD8","DC","EC","Fib","Glial","ILC","MAST","Mono/Macro","NK")

bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(cellType),nrow = 100))
names(bhatt.dist)<-cellType
names(bhatt.dist.rand)<-cellType

for (CT in cellType){
  
  cells.Inf <- colnames(Sobj)[Sobj$ParentalCluster==CT & Sobj$Class=="N"]
  cells.T <- colnames(Sobj)[Sobj$ParentalCluster==CT & Sobj$Class=="T"]
  
  tmp<-Sobj@reductions$harmony@cell.embeddings
  
  cells.Inf.pca <- tmp[cells.Inf,]
  cells.T.pca <- tmp[cells.T,]
  
  set.seed(1)
  
  for (i in 1:100) {
    
    bhatt.dist[i,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=FALSE)
    
    bhatt.dist.rand[i,CT] <- dim_dist(embed_mat_x=cells.Inf.pca,embed_mat_y=cells.T.pca,dims_use=1:30,num_cells_sample=50,distance_metric="bhatt_dist",random_sample=TRUE)
  }
}

#####################################################
#Fig2 roe

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}


