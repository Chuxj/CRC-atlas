#######################################################
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

#######################################################
#Fig2 scenic
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


