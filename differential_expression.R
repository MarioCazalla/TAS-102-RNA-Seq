library(data.table)
library(DESeq2)
library(ggplot2)
setwd("/media/mario/My Passport/IRBLLEIDA/RNA-Sequencing/alignment/")


Rmatrix <- as.matrix(read.csv("Rmatrix.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
colnames(Rmatrix) <- c("RNA-29", "RNA-30", "RNA-31", "RNA-32")

metaData <- as.matrix(read.csv("metaData.csv", header = TRUE, sep = "\t"))

dds <- DESeqDataSetFromMatrix(countData = Rmatrix,
                              colData = metaData,
                              design = ~dex, tidy = FALSE)
#Con esto nos quedamos con los genes que tienen mas de 10 counts. Disminuimos matriz
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

res <- res[order(res$padj), ]
head(res)



plotCounts(dds, gene = "CSNK1D", intgroup = "dex")
plotCounts(dds, gene = "LOC101929185", intgroup = "dex")
plotCounts(dds, gene = "GRIP2", intgroup = "dex")
plotCounts(dds, gene = "HOXB-AS3", intgroup = "dex")
plotCounts(dds, gene = "TRNS1", intgroup = "dex")
plotCounts(dds, gene = "GMPS", intgroup = "dex")


sum(res$pvalue < 0.05, na.rm = TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm = TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[order(resSig$log2FoldChange, decreasing = TRUE), ])

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
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<0.01 & abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col = "red"))

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex") #using the DESEQ2 plotPCA fxn we can

library(volcano3D)
