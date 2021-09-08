# 
# BiocManager::install(c("Biobase", "limma"))
# install.packages("devtools")
# 
# install:latest version
# devtools::install_github("Lothelab/CMScaller")
# BiocManager::install("Rtools")

library(Biobase)
library(CMScaller)
setwd("E:\\IRBLLEIDA/RNA-Sequencing/DEG_Master_Table/RECIST")
counts <- read.csv("RmatrixCDvsPD.txt", sep = "\t")
ENSEMBLids <- AnnotationDbi::select(org.Hs.eg.db,
                                              keys = counts$Geneid,
                                              keytype = "SYMBOL",
                                              columns = "ENSEMBL")
hola <- ENSEMBLids[!is.na(ENSEMBLids$ENSEMBL), ]
hola <- distinct(.data = hola, SYMBOL, .keep_all = TRUE)
adios <- filter(counts, Geneid %in% hola$SYMBOL)
adios <- cbind(adios, hola$ENSEMBL)

names <- names(adios)

par(mfrow=c(1,2))
res <- CMScaller(emat=counts, RNAseq=TRUE, FDR=0.05)


