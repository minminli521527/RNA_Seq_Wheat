### 1 Extract differentially expressed genes
### (1) DEseq2
# chrA vs chrB
allDEG2 <- read.csv('DEG_chrA_chrB_deseq2.results.csv',header = T , stringsAsFactors = F)
# padj-> P value after correction
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$log2FoldChange),]
diff_signif_DEG=na.omit(diff_signif_DEG)                       # Remove the line containing NA
# Add change column to mark genes up and down
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$padj < 0.05 & abs(diff_signif_DEG$log2FoldChange) > foldChange,
         ifelse(diff_signif_DEG$log2FoldChange > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA_chrB_deseq2_diff_signif.results.csv")

# chrA vs chrD
allDEG2 <- read.csv('DEG_chrA_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$log2FoldChange),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$padj < 0.05 & abs(diff_signif_DEG$log2FoldChange) > foldChange,
         ifelse(diff_signif_DEG$log2FoldChange > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA_chrD_deseq2_diff_signif.results.csv")

# chrB vs chrD
allDEG2 <- read.csv('DEG_chrB_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$log2FoldChange),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$padj < 0.05 & abs(diff_signif_DEG$log2FoldChange) > foldChange,
         ifelse(diff_signif_DEG$log2FoldChange > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrB_chrD_deseq2_diff_signif.results.csv")


### (2) edgeR
# chrA vs chrB
allDEG2 <- read.csv('DEG_chrA_chrB_edgeR.csv',header = T , stringsAsFactors = F)
# padj-> P value after correction
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$FDR < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)                       # Remove the line containing NA
# Add change column to mark genes up and down
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$FDR < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA_chrB_edgeR_diff_signif.csv")

# chrA vs chrD
allDEG2 <- read.csv('DEG_chrA_chrD_edgeR.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$FDR < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$FDR < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA_chrD_edgeR_diff_signif.csv")

# chrB vs chrD
allDEG2 <- read.csv('DEG_chrB_chrD_edgeR.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$FDR < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$FDR < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrB_chrD_edgeR_diff_signif.csv")


### (3) limma_voom
# chrA vs chrB
allDEG2 <- read.csv('DEG_chrA-chrB_limma_voom.csv',header = T , stringsAsFactors = F)
# padj-> P value after correction
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$adj.P.Val < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)                       # Remove the line containing NA
# Add change column to mark genes up and down
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$adj.P.Val < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA-chrB_limma_voom_diff_signif.csv")

# chrA vs chrD
allDEG2 <- read.csv('DEG_chrA-chrD_limma_voom.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$adj.P.Val < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$adj.P.Val < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrA-chrD_limma_voom_diff_signif.csv")

# chrB vs chrD
allDEG2 <- read.csv('DEG_chrB-chrD_limma_voom.csv',header = T , stringsAsFactors = F)
padj = 0.05
foldChange= 2
diff_signif_DEG= allDEG2[(allDEG2$adj.P.Val < padj & 
                          (allDEG2$logFC>foldChange | allDEG2$logFC<(-foldChange))),]
diff_signif_DEG = diff_signif_DEG[order(diff_signif_DEG$logFC),]
diff_signif_DEG=na.omit(diff_signif_DEG)
diff_signif_DEG$change = as.factor(
  ifelse(diff_signif_DEG$adj.P.Val < 0.05 & abs(diff_signif_DEG$logFC) > foldChange,
         ifelse(diff_signif_DEG$logFC > foldChange ,'UP','DOWN'),'NOT')
)
head(diff_signif_DEG)
write.csv(diff_signif_DEG,"DEG_chrB-chrD_limma_voom_diff_signif.csv")



### 2 Differentially expressed genes comparison
# Extract the list of gene names
awk -F "," '{print $2}' DEG_chrA_chrB_deseq2_diff_signif_kallisto.results.csv > AB.deseq2_kallisto.csv
awk -F "," '{print $2}' DEG_chrA_chrB_deseq2_diff_signif_last.results.csv > AB.deseq2_last.csv
awk -F "," '{print $2}' DEG_chrA_chrB_edgeR_diff_signif_kallisto.csv > AB.edgeR_kallisto.csv
awk -F "," '{print $2}' DEG_chrA_chrB_edgeR_diff_signif_last.csv > AB.edgeR_last.csv
awk -F "," '{print $2}' DEG_chrA-chrB_limma_voom_diff_signif_kallisto.csv > AB.limma_voom_kallisto.csv
awk -F "," '{print $2}' DEG_chrA-chrB_limma_voom_diff_signif_last.csv > AB.limma_voom_last.csv
awk -F "," '{print $2}' DEG_chrA_chrD_deseq2_diff_signif_kallisto.results.csv > AD.deseq2_kallisto.csv
awk -F "," '{print $2}' DEG_chrA_chrD_deseq2_diff_signif_last.results.csv > AD.deseq2_last.csv
awk -F "," '{print $2}' DEG_chrA_chrD_edgeR_diff_signif_kallisto.csv > AD.edgeR_kallisto.csv
awk -F "," '{print $2}' DEG_chrA_chrD_edgeR_diff_signif_last.csv > AD.edgeR_last.csv
awk -F "," '{print $2}' DEG_chrA-chrD_limma_voom_diff_signif_kallisto.csv > AD.limma_voom_kallisto.csv
awk -F "," '{print $2}' DEG_chrA-chrD_limma_voom_diff_signif_last.csv > AD.limma_voom_last.csv
awk -F "," '{print $2}' DEG_chrB_chrD_deseq2_diff_signif_kallisto.results.csv > BD.deseq2_kallisto.csv
awk -F "," '{print $2}' DEG_chrB_chrD_deseq2_diff_signif_last.results.csv > BD.deseq2_last.csv
awk -F "," '{print $2}' DEG_chrB_chrD_edgeR_diff_signif_kallisto.csv > BD.edgeR_kallisto.csv
awk -F "," '{print $2}' DEG_chrB_chrD_edgeR_diff_signif_last.csv > BD.edgeR_last.csv
awk -F "," '{print $2}' DEG_chrB-chrD_limma_voom_diff_signif_kallisto.csv > BD.limma_voom_kallisto.csv
awk -F "," '{print $2}' DEG_chrB-chrD_limma_voom_diff_signif_last.csv > BD.limma_voom_last.csv
awk -F "," '{print $1}' DEG_chrA-chrB_EBSeq_rsem.csv > AB.EBSeq_rsem.csv
awk -F "," '{print $1}' DEG_chrA-chrD_EBSeq_rsem.csv > AD.EBSeq_rsem.csv
awk -F "," '{print $1}' DEG_chrB-chrD_EBSeq_rsem.csv > BD.EBSeq_rsem.csv

# Install package / Load package
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(futile.logger)
library(VennDiagram)

# Differentially expressed genes vs. Venn diagram
AB.deseq2_kallisto <- read.csv('AB.deseq2_kallisto.csv',header = T , stringsAsFactors = F)
AB.deseq2_last <- read.csv('AB.deseq2_last.csv',header = T , stringsAsFactors = F)
AB.edgeR_kallisto <- read.csv('AB.edgeR_kallisto.csv',header = T , stringsAsFactors = F)
AB.edgeR_last <- read.csv('AB.edgeR_last.csv',header = T , stringsAsFactors = F)
AB.limma_voom_kallisto <- read.csv('AB.limma_voom_kallisto.csv',header = T , stringsAsFactors = F)
AB.limma_voom_last <- read.csv('AB.limma_voom_last.csv',header = T , stringsAsFactors = F)
AD.deseq2_kallisto <- read.csv('AD.deseq2_kallisto.csv',header = T , stringsAsFactors = F)
AD.deseq2_last <- read.csv('AD.deseq2_last.csv',header = T , stringsAsFactors = F)
AD.edgeR_kallisto <- read.csv('AD.edgeR_kallisto.csv',header = T , stringsAsFactors = F)
AD.edgeR_last <- read.csv('AD.edgeR_last.csv',header = T , stringsAsFactors = F)
AD.limma_voom_kallisto <- read.csv('AD.limma_voom_kallisto.csv',header = T , stringsAsFactors = F)
AD.limma_voom_last <- read.csv('AD.limma_voom_last.csv',header = T , stringsAsFactors = F)
BD.deseq2_kallisto <- read.csv('BD.deseq2_kallisto.csv',header = T , stringsAsFactors = F)
BD.deseq2_last <- read.csv('BD.deseq2_last.csv',header = T , stringsAsFactors = F)
BD.edgeR_kallisto <- read.csv('BD.edgeR_kallisto.csv',header = T , stringsAsFactors = F)
BD.edgeR_last <- read.csv('BD.edgeR_last.csv',header = T , stringsAsFactors = F)
BD.limma_voom_kallisto <- read.csv('BD.limma_voom_kallisto.csv',header = T , stringsAsFactors = F)
BD.limma_voom_last <- read.csv('BD.limma_voom_last.csv',header = T , stringsAsFactors = F)

AB.EBSeq_rsem <- read.csv('AB.EBSeq_rsem.csv',header = T , stringsAsFactors = F)
AD.EBSeq_rsem <- read.csv('AD.EBSeq_rsem.csv',header = T , stringsAsFactors = F)
BD.EBSeq_rsem <- read.csv('BD.EBSeq_rsem.csv',header = T , stringsAsFactors = F)

# beautify fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),alpha = 0.4,cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4")
# col = "black", Border color
# fontface = "bold", label font
# fill = c("cornflowerblue", "green", "yellow", "darkorchid1"), fill color
# alpha = 0.8, transparency
# cex = 4, label font size
# cat.cex = 3, the size of the class name
# cat.fontface = "bold", class name body
# margin = 0.04, marginal distance


# 2
venn.diagram(list(AB.deseq2_kallisto=AB.deseq2_kallisto$X,AB.deseq2_last=AB.deseq2_last$X),"AB.deseq2_kallisto=deseq2_last.pdf",col = "black",fill = c("cornflowerblue", "green"),alpha = 0.8,cat.col = c("darkblue", "darkgreen"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)

# 2
venn.diagram(list(AB.deseq2_kallisto=AB.deseq2_kallisto$X,AB.EBSeq_rsem=AB.EBSeq_rsem$Pattern1),"AB.deseq2_kallisto=EBSeq_rsem.pdf",col = "black",fill = c("cornflowerblue", "green"),alpha = 0.8,cat.col = c("darkblue", "darkgreen"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)

# 4
venn.diagram(list(AB.deseq2_kallisto=AB.deseq2_kallisto$X,AB.deseq2_last=AB.deseq2_last$X,AB.edgeR_kallisto=AB.edgeR_kallisto$X, AB.edgeR_last=AB.edgeR_last$X),"AB.deseq2_kallisto=deseq2_last=edgeR_kallisto=edgeR_last.pdf",col = "black",fill = c("cornflowerblue", "green","yellow", "darkorchid1"),alpha = 0.8,cat.col = c("darkblue", "darkgreen","orange", "darkorchid4"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)

# 4
venn.diagram(list(AD.deseq2_kallisto=AD.deseq2_kallisto$X,AD.deseq2_last=AD.deseq2_last$X,AD.edgeR_kallisto=AD.edgeR_kallisto$X, AD.edgeR_last=AD.edgeR_last$X),"AD.deseq2_kallisto=deseq2_last=edgeR_kallisto=edgeR_last.pdf",col = "black",fill = c("cornflowerblue", "green","yellow", "darkorchid1"),alpha = 0.8,cat.col = c("darkblue", "darkgreen","orange", "darkorchid4"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)

# 4
venn.diagram(list(BD.deseq2_kallisto=BD.deseq2_kallisto$X,BD.deseq2_last=BD.deseq2_last$X,BD.edgeR_kallisto=BD.edgeR_kallisto$X, BD.edgeR_last=BD.edgeR_last$X),"BD.deseq2_kallisto=deseq2_last=edgeR_kallisto=edgeR_last.pdf",col = "black",fill = c("cornflowerblue", "green","yellow", "darkorchid1"),alpha = 0.8,cat.col = c("darkblue", "darkgreen","orange", "darkorchid4"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)

# 5
venn.diagram(list(AB.deseq2_kallisto=AB.deseq2_kallisto$X,AB.deseq2_last=AB.deseq2_last$X,AB.EBSeq_rsem=AB.EBSeq_rsem$Pattern1,AB.edgeR_kallisto=AB.edgeR_kallisto$X, AB.edgeR_last=AB.edgeR_last$X),"AB.deseq2_kallisto=deseq2_last=EBSeq_rsem=edgeR_kallisto=edgeR_last.pdf",col = "black",fill = c("cornflowerblue", "green","yellow", "darkorchid1","red"),alpha = 0.8,cat.col = c("darkblue", "darkgreen","orange", "darkorchid4","darkred"), cat.cex = 1.0,rotation.degree = 0, fontface = "bold",margin = 0.2)


