### dependencies: run if not already installed
### limma has lof of dependencies - takes time
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "limma"))
# install.packages("devtools")

### install: latest version
devtools::install_github("Lothelab/CMScaller")
install.packages("Rtools")
library(Biobase)
library(CMScaller)
par(mfrow=c(1,2))

### CMS prediction of TCGA primary colorectal cancers
counts <- exprs(crcTCGAsubset)
head(counts)
class(crcTCGAsubset)
class(counts)
res <- CMScaller(crcTCGAsubset, RNAseq=TRUE, doPlot=TRUE)
head(res)

### Camera Gene Set Analysis with CMS informative gene sets
cam <- CMSgsa(emat=crcTCGAsubset, class=res$prediction, RNAseq=TRUE)
head(cam$CMS4)

### limma differential gene expression analysis and visualization
deg <- subDEG(emat=crcTCGAsubset, class=res$prediction, doVoom=TRUE)
subVolcano(deg, geneID="symbol")
