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

setwd("C:\\Users/usuari/Desktop/TAS-102-RNA-Seq/CMScaller")

counts <- read.csv("RmatrixCMScaller.txt", sep = "\t")

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
colnames(adios)[90] <- "entrez"
names <- names(adios)
names

encantado <- adios[ , c(90,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78)]
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
write.table(res, "CMS_results.csv", sep = "\t")
head(res)
hist(res$p.value)

cam <- CMSgsa(emat=x, class=res$prediction, RNAseq=TRUE)
head(cam)