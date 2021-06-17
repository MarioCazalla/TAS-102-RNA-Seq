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

setwd("/media/mario/My Passport/IRBLLEIDA/RNA-S")
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/")

#Preparamos la matriz de analisis para DESeq2

Rmatrix <- as.matrix(read.csv("Rmatrix.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-29", "RNA-30", "RNA-31", "RNA-32")
dim(Rmatrix)
metaData <- as.matrix(read.csv("metaData.csv", header = TRUE, sep = "\t"))

dds <- DESeqDataSetFromMatrix(countData = Rmatrix,
                              colData = metaData,
                              design = ~dex, tidy = FALSE)
#Con esto nos quedamos con los genes que tienen mas de 1 counts. Disminuimos matriz
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
# nrow(dds)

#DESeq2
dds <- DESeq(dds)
res <- results(dds)

#Observamos los resultados
head(results(dds, tidy=TRUE))
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
#Vemos los resultados
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = FALSE), ])

#Before multiple testing: genes with pvalue < 0.05
resSIG <- subset(res, pvalue < 0.05)
#Vemos los resultados
head(resSIG[order(resSIG$log2FoldChange, decreasing = TRUE), ])

#Ahora filtraremos por log2FoldChange para ver cuales son UP o DOWN regulated
  #ADJUSTED P-VALUE
log2FC_UP <- subset(resSig, resSig$log2FoldChange > 0) #3
log2FC_UP <- log2FC_UP[order(log2FC_UP$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN <- subset(resSig, resSig$log2FoldChange < 0) #3
log2FC_DOWN <- log2FC_DOWN[order(log2FC_DOWN$log2FoldChange, decreasing = FALSE), ]

resSig_up_names <- log2FC_UP@rownames
resSig_down_names <- log2FC_DOWN@rownames
  #P-VALUE
log2FC_UP_pv <- subset(resSIG, resSIG$log2FoldChange > 0) #41
log2FC_UP_pv <- log2FC_UP_pv[order(log2FC_UP_pv$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN_pv <- subset(resSIG, resSIG$log2FoldChange < 0) #24
log2FC_DOWN_pv <- log2FC_DOWN_pv[order(log2FC_DOWN_pv$log2FoldChange, decreasing = FALSE), ]

resSIG_UP_names <-log2FC_UP_pv@rownames
resSIG_DOWN_names <- log2FC_DOWN_pv@rownames
#Guardamos los resultados 
write.table(resSig_down_names, file = "down_regulated_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_up_names, file = "up_regulated_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(resSIG_DOWN_names, file = "down_regulated_pvalue_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_UP_names, file = "up_regulated_pvalue_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#Con esto vamos a generar la lista de KEGG pathways
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = resSIG_DOWN_names, #Cambiar key segun analisis
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")

entrez_ids <- entrez_ids[ , "ENTREZID"] 
entrez_ids <- entrez_ids[!is.na(entrez_ids)] 
KEGG_analysis <-enrichKEGG(entrez_ids,
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


#Creamos los plots en los que vemos las diferencias entre Resistente y parental
resSIG[order(resSIG$log2FoldChange, decreasing = TRUE), ]

geneCounts <- plotCounts(dds, gene = "HOXB-AS3", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "LOC101929185", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GRIP2", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "HOXB-AS3", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GMPS", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)


geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = dex, y = count, color = celltype, group = celltype)) + scale_y_log10() + geom_point(size=3) + geom_line()


#plotCounts del gen con mayos p.adj
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = ("dex"))


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resSig, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resSig, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resSig, padj<0.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col = "red"))

