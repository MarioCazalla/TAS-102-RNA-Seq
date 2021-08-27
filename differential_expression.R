## Install BiocManager if not installed ##

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.13")
## Install packages if not installed ##
#BiocManager::install(c("EnhancedVolcano", "ReactomePA", "data.table", "dplyr", "DESeq2", "ggplot2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "ggbeeswarm", "volcano3D"))

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
library(EnhancedVolcano)

#library(volcano3D)

##### SNU-C4 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/SNU-C4/")
##### LS513 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/LS513/")
##### LIM2099 #####
setwd("E:\\IRBLLEIDA/RNA-Sequencing/LIM-2099/")

### Preparamos la matriz de analisis para DESeq2 ###

##### SNU-C4 #####
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

### 6 vs 6 MasterTable ### Analisis de las 6 lineas celulares de cada extremo.
setwd("C:\\Users/usuari/Desktop/TAS-102-RNA-Seq/DEG_Master_Table/6vs6/")
Rmatrix <- as.matrix(read.csv("Rmatrix6vs6.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
metaData <- as.matrix(read.csv("metaData6vs6.txt", header = TRUE, sep = "\t"))

### 6 vs Rest MasterTable ### Analysis of 6 cell lines of the right side of histogram vs ALL ###
setwd("C:\\Users/usuari/Desktop/TAS-102-RNA-Seq/DEG_Master_Table/6vsALL/")
Rmatrix <- as.matrix(read.csv("Rmatrix6vsALL.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
metaData <- as.matrix(read.csv("metaData6vsALL.txt", header = TRUE, sep = "\t"))

### CD vs PD ### Analysis of 6 cell lines of the right side of histogram vs ALL ###
setwd("C:\\Users/usuari/Desktop/TAS-102-RNA-Seq/DEG_Master_Table/RECIST/")
Rmatrix <- as.matrix(read.csv("RmatrixCDvsPD.txt", header = TRUE, sep = "\t", row.names = "Geneid"))
metaData <- as.matrix(read.csv("metaDataCDvsPD.txt", header = TRUE, sep = "\t"))

###############

Rmatrix1 <- na.omit(Rmatrix)
dim(Rmatrix1)
dds <- DESeqDataSetFromMatrix(countData = Rmatrix1,
                              colData = metaData,
                              design = ~dex, tidy = FALSE)

### This command let us to save genes with more than 1 count, matrix is smaller ###
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

### Run DESeq2 ###
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

### Sort by p.adj ###
res <- res[order(res$padj), ]
head(res)

### Filtering: multiple test correction: genes with padj < 0.2 ###
resSig <- subset(res, padj < 0.2)
head(resSig[ order(resSig$log2FoldChange, decreasing = FALSE), ]) #To see results

#### Calcular el porcentaje de genes que han cambiado tras el tratamiento



### Save dataframe wit gene name, log2FC, pvalue and padj
resSig_log2FC <- as.data.frame(cbind(resSig@rownames, resSig$log2FoldChange, resSig$pvalue, resSig$padj))
colnames(resSig_log2FC) <- c("GeneName", "log2FC", "pvalue", "padj")
write.table(resSig_log2FC, file = "genename_log2FC.csv", sep = "\t", row.names = FALSE) #Next step: convert "." to "," in excel to treat them like numbers, not like text.

### Now, we are filtering by log2FoldChange to clasify which are up or down regulated.
  #Using the data that pass multiple test correction (ADJUSTED P-VALUE < 0.02)
log2FC_UP <- subset(resSig, resSig$log2FoldChange > 1.5) 
log2FC_UP <- log2FC_UP[order(log2FC_UP$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN <- subset(resSig, resSig$log2FoldChange < -1.5) 
log2FC_DOWN <- log2FC_DOWN[order(log2FC_DOWN$log2FoldChange, decreasing = FALSE), ]

#We save the downregulated gene names with their log2fc.
name_log2fc_down <- as.data.frame(log2FC_DOWN@rownames)
colnames(name_log2fc_down) <- "GeneName"
name_log2fc_down$log2fc=log2FC_DOWN$log2FoldChange
write.table(name_log2fc_down, file = "name_log2FC_down.csv", sep = "\t", row.names = FALSE)

#We save the upregulated gene names with their log2fc.
name_log2fc_up <- as.data.frame(log2FC_UP@rownames)
colnames(name_log2fc_up) <- "GeneName"
name_log2fc_up$log2fc=log2FC_UP$log2FoldChange
write.table(name_log2fc_up, file = "name_log2FC_up.csv", sep = "\t", row.names = FALSE)

### We keep gene names
resSig_names <- resSig@rownames
resSig_up_names <- log2FC_UP@rownames
resSig_down_names <- log2FC_DOWN@rownames

### Before multiple test correction: genes with pvalue < 0.05
resSIG <- subset(res, pvalue < 0.05)
head(resSIG[order(resSIG$log2FoldChange, decreasing = TRUE), ])

  #Using the data that pass P-VALUE < 0.05
log2FC_UP_pv <- subset(resSIG, resSIG$log2FoldChange > 0) 
log2FC_UP_pv <- log2FC_UP_pv[order(log2FC_UP_pv$log2FoldChange, decreasing = TRUE), ]
log2FC_DOWN_pv <- subset(resSIG, resSIG$log2FoldChange < 0) 
log2FC_DOWN_pv <- log2FC_DOWN_pv[order(log2FC_DOWN_pv$log2FoldChange, decreasing = FALSE), ]

#We save the downregulated gene names with their log2fc.
name_log2fc_down_pv <- as.data.frame(log2FC_DOWN_pv@rownames)
colnames(name_log2fc_down_pv) <- "GeneName"
name_log2fc_down_pv$log2fc=log2FC_DOWN_pv$log2FoldChange
write.table(name_log2fc_down_pv, file = "name_log2FC_down.csv", sep = "\t", row.names = FALSE)

#We save the upregulated gene names with their log2fc.
name_log2fc_up_pv <- as.data.frame(log2FC_UP_pv@rownames)
colnames(name_log2fc_up_pv) <- "GeneName"
name_log2fc_up_pv$log2fc=log2FC_UP_pv$log2FoldChange
write.table(name_log2fc_up_pv, file = "name_log2FC_up.csv", sep = "\t", row.names = FALSE)

### We keep gene names
resSIG_names <- resSIG@rownames
resSIG_UP_names <-log2FC_UP_pv@rownames
resSIG_DOWN_names <- log2FC_DOWN_pv@rownames

### To remove the gene names LOC... ### 
#no_LOC <- resSIG[!grepl("^LOC", resSIG@rownames),]

### To keep the results in the directory ###
write.table(resSig_names, file = "diff_expressed_padj_genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_down_names, file = "down_regulated_padj_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSig_up_names, file = "up_regulated_padj_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

write.table(resSIG_names, file = "diff_expressed_pval__genenames.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_DOWN_names, file = "down_regulated_pval_names.csv", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(resSIG_UP_names, file = "up_regulated_pval_names.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#### UNIPROT ID CONVERSION FROM x ####

UNIPROT_IDs <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = resSig@rownames, #Change Keys when you change the analysis
                                     keytype = "SYMBOL",
                                     columns = "ENSEMBL") #Change Keys when you change the analysis
UNIPROT_IDs <- UNIPROT_IDs[!is.na(UNIPROT_IDs$UNIPROT), ]
write.table(UNIPROT_IDs, file = "uniprots_id.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#### CYTOSCAPE dataframe Preparation ####
cytoscape <- as.data.frame(resSig@rownames)
colnames(cytoscape) <- "Gene_Name"
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = resSig@rownames, #Change Keys when you change the analysis
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")
cytoscape$GeneID=entrez_ids$ENTREZID
cytoscape$log2FC=resSig$log2FoldChange
cytoscape$pvalue=resSig$pvalue
cytoscape$Gene_Name <- NULL
cytoscape <- cytoscape[!is.na(cytoscape$GeneID), ] 
write.table(cytoscape, file = "cytoscape_ENTREZID.csv", sep = "\t", row.names = FALSE)

#### KEGG pathways ####
entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = resSig_names, #Change Keys when you change the analysis
                                    keytype = "SYMBOL",
                                    columns = "ENTREZID")

## We create a dataframe with EntrezID and GeneName to help us to know the genename after doing enrichKEGG##
name_entrezid <- entrez_ids

entrez_ids <- entrez_ids[ , "ENTREZID"] 
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
### Biological Process ###
GO_BP <- GO_analysis(resSIG_names, "BP")
GO_BP_up <- GO_analysis(resSig_up_names, "BP")
GO_BP_down <- GO_analysis(resSig_down_names, "BP")

## We save the results in a dataframe to be able to investigate them ###
GO_BP_df <- as.data.frame(GO_BP@result)
GO_BP_up_df <- as.data.frame(GO_BP_up@result)
GO_BP_down_df <- as.data.frame(GO_BP_down@result)

### Molecular Function and Cellular components ###
GO_MF <- GO_analysis(resSig_names, "MF")
GO_CC <- GO_analysis(resSig_names, "CC")

## We save the results in a dataframe to be able to investigate them ###
GO_MF_df <- as.data.frame(GO_MF@result)
GO_CC_df <- as.data.frame(GO_CC@result)

### We save the results in a .csv file ###
write.table(GO_BP_df, file = "GO_BP.csv", sep = "\t", row.names = FALSE)
write.table(GO_BP_up_df, file = "GO_BP_up.csv", sep = "\t", row.names = FALSE)
write.table(GO_BP_down_df, file = "GO_BP_down.csv", sep = "\t", row.names = FALSE)
write.table(GO_MF_df, file = "GO_MF_pval.csv", sep = "\t")
write.table(GO_CC_df, file = "GO_CC_pval.csv", sep = "\t")

#We create the plots in which we can see the differences between Resistants and Sensibles ###
geneCounts <- plotCounts(dds, gene = "TM4SF4", intgroup = c("dex", "celltype", "id"), returnData=TRUE)

# Plotting the X gene normalized counts, using the samplenames (rownames of d as labels)
ggplot(geneCounts, aes(x = dex, y = count, color = celltype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(geneCounts))) + 
  theme_bw() +
  ggtitle("TM4SF4") +
  theme(plot.title = element_text(hjust = 0.5))

### To see the expression of a gene in each cell lines ###
d <- plotCounts(dds, gene="TM4SF18", intgroup="dex", returnData=TRUE)
d$name <- rownames(d)
ggplot(d, aes(x=dex, y=count, color=dex)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  geom_text_repel(aes(label = name)) + 
  theme_bw() +
  ggtitle("Gene KRT7") +
  theme(plot.title=element_text(hjust=0.5))

### VOLCANO PLOT to see easier the log2fc and the p value ###
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', # change if we use padj
                title = 'TAS102 6 vs ALL',
                pCutoff = 0.05,
                FCcutoff = 2.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
