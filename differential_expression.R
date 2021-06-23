# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.13")
# 
# BiocManager::install(c("data.table", "dplyr", "DESeq2", "ggplot2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "ggbeeswarm", "volcano3D"))
#install.packages("countToFPKM")
library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggbeeswarm)
library(volcano3D)

#Linux
setwd("/media/mario/My Passport/IRBLLEIDA/RNA-Sequencing/alignment/SNU-C4")
#Windows
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/LIM-2099//")

#Preparamos la matriz de analisis para DESeq2

Rmatrix <- as.matrix(read.csv("Rmatrix_lim2099.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-25", "RNA-26", "RNA-27", "RNA-28")
metaData <- as.matrix(read.csv("metaData_lim2099.csv", header = TRUE, sep = "\t"))

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
resSIG <- subset(res, pvalue < 0.06)
#Vemos los resultados
head(resSIG[order(resSIG$log2FoldChange, decreasing = TRUE), ])

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

#Guardamos los resultados 
write.table(resSig_names, file = "diff_expressed_padj_genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_down_names, file = "down_regulated_padj_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_up_names, file = "up_regulated_padj_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(resSIG_names, file = "diff_expressed_pval__genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_DOWN_names, file = "down_regulated_pval_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_UP_names, file = "up_regulated_pval_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

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

write.table(GO_BP_df, file = "GO_BP_padj.csv", sep = "\t")
write.table(GO_MF_df, file = "GO_MF_padj.csv", sep = "\t")
write.table(GO_CC_df, file = "GO_CC_padj.csv", sep = "\t")


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


