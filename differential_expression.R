if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

BiocManager::install(c("DESeq2", "ggplot2"))
library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
setwd("/media/mario/My Passport/IRBLLEIDA/RNA-S")
setwd("E:\\IRBLLEIDA/RNA-Sequencing/alignment/")


Rmatrix <- as.matrix(read.csv("Rmatrix.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-29", "RNA-30", "RNA-31", "RNA-32")

metaData <- as.matrix(read.csv("metaData.csv", header = TRUE, sep = "\t"))

dds <- DESeqDataSetFromMatrix(countData = Rmatrix,
                              colData = metaData,
                              design = ~dex, tidy = FALSE)
#Con esto nos quedamos con los genes que tienen mas de 10 counts. Disminuimos matriz
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

#install.packages("countToFPKM")

res <- res[order(res$padj), ]
head(res)

sum(res$pvalue < 0.05, na.rm = TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm = TRUE)
#After multiple testing 
resSig <- subset(res, padj < 0.1)
#Before multiple testing
resSIG <- subset(res, pvalue < 0.05)
head(resSig[order(resSig$log2FoldChange, decreasing = TRUE), ])
#ADJUSTED P-VALUE
log2FC_UP <- subset(resSig, resSig$log2FoldChange > 1) #3
log2FC_UP@rownames
log2FC_DOWN <- subset(resSig, resSig$log2FoldChange < -1) #3
log2FC_DOWN@rownames
#P-VALUE
log2FC_UP_pv <- subset(resSIG, resSIG$log2FoldChange > 1) #41
log2FC_UP_pv@rownames
log2FC_DOWN_pv <- subset(resSIG, resSIG$log2FoldChange < -1) #24
log2FC_DOWN_pv@rownames

resSig_names <- resSig@rownames
resSIG_names <- resSIG@rownames

library(ggbeeswarm)
geneCounts <- plotCounts(dds, gene = "CSNK1D", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "LOC101929185", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GRIP2", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "HOXB-AS3", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

geneCounts <- plotCounts(dds, gene = "GMPS", intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)


#plotCounts
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = ("dex"))

library(ggbeeswarm)
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex", "celltype"), returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + scale_y_log10() + geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = dex, y = count, color = celltype, group = celltype)) + scale_y_log10() + geom_point(size=3) + geom_line()


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resSig, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resSig, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resSig, padj<0.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col = "red"))

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex") #using the DESEQ2 plotPCA fxn we can

library(volcano3D)
