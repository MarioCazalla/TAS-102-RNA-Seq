# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.13")
# 
# BiocManager::install(c("ReactomePA", "data.table", "dplyr", "DESeq2", "ggplot2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "ggbeeswarm", "volcano3D"))
#install.packages("countToFPKM")
library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggbeeswarm)
library(DOSE)
library(ReactomePA)
#library(volcano3D)

######  SNU-C4 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/SNU-C4/")
##### LS513 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/LS513//")
##### LIM2099 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/LIM-2099//")

#Preparamos la matriz de analisis para DESeq2

######  SNU-C4 #####
Rmatrix <- as.matrix(read.csv("Rmatrix_SNUC4.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-29", "RNA-30", "RNA-31", "RNA-32") 
metaData <- as.matrix(read.csv("metaData_SNUC4.csv", header = TRUE, sep = "\t"))

##### LS513 #####
Rmatrix <- as.matrix(read.csv("Rmatrix_LS513.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-33", "RNA-34", "RNA-35", "RNA-36") 
metaData <- as.matrix(read.csv("metaData_LS513.csv", header = TRUE, sep = "\t"))

##### LIM2099 #####
Rmatrix <- as.matrix(read.csv("Rmatrix_LIM2099.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-25", "RNA-26", "RNA-27", "RNA-28") 
metaData <- as.matrix(read.csv("metaData_LIM2099.csv", header = TRUE, sep = "\t"))


dim(Rmatrix)

Rmatrix1 <- na.omit(Rmatrix)
dim(Rmatrix1)
dds <- DESeqDataSetFromMatrix(countData = Rmatrix1,
                              colData = metaData,
                              design = ~dex, tidy = FALSE)
#Con esto nos quedamos con los genes que tienen mas de 1 counts. Disminuimos matriz
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
# nrow(dds)

#DESeq2

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

#rm(keep, metaData, Rmatrix, dds)
#Observamos los resultados

summary(res)
#Ordenamos por p.adj
res <- res[order(res$padj), ]
#Vemos los resultados
head(res)

# sum(res$pvalue < 0.05, na.rm = TRUE)
# sum(!is.na(res$pvalue))
# sum(res$padj < 0.1, na.rm = TRUE)

#After multiple testing: genes with padj < 0.2 
resSig <- subset(res, padj < 0.2)

#### Calcular el porcentaje de genes que han cambiado tras el tratamiento



#Vemos los resultados
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = FALSE), ])

#Before multiple testing: genes with pvalue < 0.05
resSIG <- subset(res, pvalue < 0.05)
#Vemos los resultados
head(resSIG[order(resSIG$log2FoldChange, decreasing = TRUE), ])

#Guardamos gene name con log2FC y pvalue/padj

resSIG_log2FC <- as.data.frame(cbind(resSIG@rownames, resSIG$log2FoldChange, resSIG$pvalue, resSIG$padj))
colnames(resSIG_log2FC) <- c("GeneName", "log2FC", "pvalue", "padj")
write.table(resSIG_log2FC, file = "genename_log2FC.csv", sep = "\t", row.names = FALSE)

#Ahora filtraremos por log2FoldChange para ver cuales son UP o DOWN regulated
  #ADJUSTED P-VALUE
log2FC_UP <- subset(resSig, resSig$log2FoldChange > 0) 
log2FC_UP <- log2FC_UP[order(log2FC_UP$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN <- subset(resSig, resSig$log2FoldChange < 0) 
log2FC_DOWN <- log2FC_DOWN[order(log2FC_DOWN$log2FoldChange, decreasing = FALSE), ]

resSig_names <- resSig@rownames
resSig_up_names <- log2FC_UP@rownames
resSig_down_names <- log2FC_DOWN@rownames

  #P-VALUE
log2FC_UP_pv <- subset(resSIG, resSIG$log2FoldChange > 0) 
log2FC_UP_pv <- log2FC_UP_pv[order(log2FC_UP_pv$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN_pv <- subset(resSIG, resSIG$log2FoldChange < 0) 
log2FC_DOWN_pv <- log2FC_DOWN_pv[order(log2FC_DOWN_pv$log2FoldChange, decreasing = FALSE), ]

resSIG_names <- resSIG@rownames
resSIG_UP_names <-log2FC_UP_pv@rownames
resSIG_DOWN_names <- log2FC_DOWN_pv@rownames


#### UNIPROT ####
UNIPROT_IDs <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = resSIG_names, #Cambiar key segun analisis
                                     keytype = "SYMBOL",
                                     columns = "UNIPROT")
UNIPROT_IDs <- UNIPROT_IDs[!is.na(UNIPROT_IDs$UNIPROT), ]
UNIPROT_IDs <- UNIPROT_IDs[, "UNIPROT"]
write.table(UNIPROT_IDs, file = "uniprots_id.csv", sep = ",", row.names = FALSE, col.names = FALSE)



#Guardamos los resultados 
write.table(resSig_names, file = "diff_expressed_padj_genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_down_names, file = "down_regulated_padj_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_up_names, file = "up_regulated_padj_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(resSIG_names, file = "diff_expressed_pval__genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_DOWN_names, file = "down_regulated_pval_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_UP_names, file = "up_regulated_pval_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#### CYTOSCAPE ####
cytoscape <- as.data.frame(resSIG@rownames)
colnames(cytoscape) <- "Gene_Name"
cytoscape$log2FC=resSIG$log2FoldChange
cytoscape$pvalue=resSIG$pvalue
cytoscape$GeneID=entrez_ids$ENSEMBL

### REACTOMEPA ###
no_LOC <- resSIG[!grepl("^LOC", resSIG@rownames),]

eg = bitr(no_LOC@rownames, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- eg[, "ENTREZID"]
x <- enrichPathway(gene=eg,pvalueCutoff=0.05, readable=T)
head(x)


x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)




#Para transformar a ENSG...(keytype es lo que tú quieres cambiar)(column a lo que quieres cambiar) CYTOSCAPE:
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = resSIG@rownames, #Cambiar key segun analisis
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")


entrez_ids <- entrez_ids[!is.na(entrez_ids$ENSEMBL), ]
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

write.table(cytoscape, file = "cytoscape.csv", sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(entrez_ids, file = "entrez_cytoscape.csv", sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(resSIG$pvalue, file = "pvalue_cytoscape.csv", sep = "\t", col.names = FALSE, row.names = FALSE)

#### KEGG pathways ####
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = resSIG_names, #Cambiar key segun analisis
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")

entrez_ids <- entrez_ids[ , "ENTREZID"] 
#Voy a crear dataframe con los entrezid y los nombres de los genes
name_entrezid <- entrez_ids

entrez_ids <- entrez_ids[!is.na(entrez_ids)] 
KEGG_analysis <- enrichKEGG(entrez_ids,
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           #universe,
                           minGSSize = 10,
                           maxGSSize = 500,
                           qvalueCutoff = 0.2,
                           use_internal_data = FALSE)

KEGG_pathways <- as.data.frame(KEGG_analysis@result[["ID"]])
KEGG_pathways$`KEGG_analysis@result[["ID"]]`<- NULL
KEGG_pathways$ID=KEGG_analysis@result[["ID"]]
KEGG_pathways$Description=KEGG_analysis@result[["Description"]]
KEGG_pathways$geneID=KEGG_analysis@result[["geneID"]]
KEGG_pathways$Count=KEGG_analysis@result[["Count"]]
KEGG_pathways$p_adjust=KEGG_analysis@result[["p.adjust"]]
KEGG_pathways$pvalue=KEGG_analysis@result[["pvalue"]]
KEGG_pathways$GeneRatio=KEGG_analysis@result[["GeneRatio"]]

write.table(KEGG_pathways, file = "KEGG_pathways_padj.csv", row.names = FALSE)
write.table(KEGG_pathways, file = "KEGG_pathways_pval.csv", row.names = FALSE)

#### GO:BP ####
GO_analysis <- function (genes, ontology){
  clusterProfiler::enrichGO(gene          = genes,
                            # universe      = universe,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'SYMBOL',
                            ont           = ontology, # "BP", "MF", "CC" o "ALL"
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
}

GO_BP <- GO_analysis(resSIG_names, "BP")
GO_MF <- GO_analysis(resSIG_names, "MF")
GO_CC <- GO_analysis(resSIG_names, "CC")

GO_BP_df <- as.data.frame(GO_BP@result)
GO_MF_df <- as.data.frame(GO_MF@result)
GO_CC_df <- as.data.frame(GO_CC@result)

write.table(GO_BP_df, file = "GO_BP_pval.csv", sep = "\t")
write.table(GO_MF_df, file = "GO_MF_pval.csv", sep = "\t")
write.table(GO_CC_df, file = "GO_CC_pval.csv", sep = "\t")

#Creamos los plots en los que vemos las diferencias entre Resistente y parental

  #Aqui si queremos ver diferentes colores entre lineas celulares, en color ponemos
    #celltype, si es la misma linea celular, ponemos dex para diferenciar R de Parental

geneCounts <- plotCounts(dds, gene = "MED22", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = dex)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "ZACN", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = dex)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GRIP2", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "HOXB-AS3", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GMPS", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)


geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = dex, y = count, color = celltype, group = celltype)) + scale_y_log10() + geom_point(size=3) + geom_line()


