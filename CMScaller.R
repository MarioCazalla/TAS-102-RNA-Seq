# install.packages("devtools")
# devtools::install_github("Lothelab/CMScaller")
# BiocManager::install("snow", "edgeR", "randomForest", "Biobase", "limma", "CMScaller", "dplyr", "org.Hs.eg.db")
# BiocManager::install("Rtools", force = TRUE)


library(Biobase)
library(CMScaller)
library(dplyr)
library(org.Hs.eg.db)
library(randomForest)
library(snow)
library(edgeR)

setwd("E:\\IRBLLEIDA/RNA-Sequencing/DEG_Master_Table/RECIST")

counts <- read.csv("RmatrixCDvsPD.txt", sep = "\t")

names(counts)

ENSEMBLids <- AnnotationDbi::select(org.Hs.eg.db,
                                            keys = counts$Geneid,
                                            keytype = "SYMBOL",
                                            columns = "ENTREZID")


columns(org.Hs.eg.db)
incomplete <- ENSEMBLids[!is.na(ENSEMBLids$ENTREZ), ]
incomplete <- distinct(.data = incomplete, SYMBOL, .keep_all = TRUE)
adios <- filter(counts, Geneid %in% incomplete$SYMBOL)
adios <- cbind(adios, incomplete$ENTREZID)
colnames(adios)[79] <- "entrez"
names <- names(adios)
names

encantado <- adios[ , c(79,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78)]
df <- encantado[,-1]

rownames(df) <- encantado[,1]


class(df)
df <- data.matrix(df)

head(df)
dim(df)

expresion <- ExpressionSet(assayData = df)
x <- exprs(expresion)
class(x)

par(mfrow=c(1,2))
res <- CMScaller(emat=x,rowNames = "entrez", RNAseq=TRUE, FDR=0.05)

CMScaller <- function(emat, templates=CMScaller::templates.CMS,
                      rowNames="ensg",
                      RNAseq=TRUE, nPerm=1000, seed=NULL,
                      FDR=0.05, doPlot=TRUE, verbose=TRUE) {
  
  # checkInput ##############################################################
  
  # check datatype input and try to coerce to matrix
  if (class(emat)[1] == "ExpressionSet") {
    emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
  }
  if (class(emat)[1] == "data.frame") emat <- as.matrix(emat)
  if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
  if (is.null(rownames(emat))) stop("missing Ensembl id rownames(emat)")
  
  if (ncol(emat) < 30) warnings("few samples - high prediction variance",
                                call.=FALSE)
  
  if (rowNames != "entrez") {
    if (!rowNames %in% c("symbol", "ensg"))
      stop("invalid rowNames, must be either entrez, symbol or ensg")
    emat <- replaceGeneId(emat, id.in=rowNames, id.out="entrez")
  }
  
  # log2-transform and quantile normalize RNA-seq data
  if (isTRUE(RNAseq)) {
    if (isTRUE(verbose))
      message("performing log2-transform and quantile normalization...")
    emat <- limma::normalizeQuantiles(log2(emat+.25))
  }
  
  # sanity check - whether rownames appear to be Entrez ids
  is.na.rows <- is.na(fromTo(rownames(emat), rough=TRUE))
  mm <- sum(is.na.rows)/nrow(emat)
  if (mm > 0.15) {
    message (paste0(sum(is.na.rows),"/",nrow(emat),
                    " rownames(emat) failed to match to human gene identifiers"))
    warning (paste0("verify that rownames(emat) are ", rowNames),
             call.=FALSE)
  }
  
  # scale and center data, basically a wrapper for scale() function
  emat <- ematAdjust(emat)
  
  # ntpPredict ##############################################################
  
  res <- ntp(emat, templates, seed=seed, nPerm=nPerm,
             doPlot=doPlot, verbose=verbose)
  res <- subSetNA(res, FDR=FDR, verbose=verbose)
  
  # output ##################################################################
  
  # sanity check III - whether any FDR-values are above .1
  if (nPerm > 500) if (min(res$FDR) > .1)
    warning("low-confidence predictions - check input",call.=FALSE)
  
  return(res)
}
CMScaller(x)
