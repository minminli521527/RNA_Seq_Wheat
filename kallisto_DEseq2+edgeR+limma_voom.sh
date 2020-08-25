### 1 Kallisto Quantitative
# 1.1 Index the transcript
# Option: -i indicates for index, followed by index name. If you still use this species in the future, you can directly use the index created this time without repeating this step.
# Triticum
kallisto index /kallisto/0_data/Triticum_aestivum.IWGSC.cdna.all.fa.gz -i /kallisto/0_data/Triticum_aestivum.index


# 1.2 Identification and quantification of transcripts
# Paired-end sequencing
# -b: The number of samples for the self-service method, the author recommends to be set to at least 30, if the downstream does not need to use the booststraps value (for example, do not do sleuth, QTL analysis), -b can be set to 0 to speed up the operation.

# Triticum
mkdir 1_result
cd 1_result
# Create a file soft link
ln -s '/software/0_data/Triticum_AC1-rep1_1.fastq' '/software/0_data/Triticum_AC1-rep1_2.fastq' '/software/0_data/Triticum_AC1-rep2_1.fastq' '/software/0_data/Triticum_AC1-rep2_2.fastq' '/software/0_data/Triticum_AC4-rep1_1.fastq' '/software/0_data/Triticum_AC4-rep1_2.fastq' '/software/0_data/Triticum_AC4-rep2_1.fastq' '/software/0_data/Triticum_AC4-rep2_2.fastq' ./
# Quantitative
ls Triticum*_1.fastq |while read id
do
kallisto quant --pseudobam -i ./Triticum_aestivum.index \
-o ./${id%_*} -t 4 -b 30 \
./$id ./${id%_*}_2.fastq \
| samtools view -F 4 -Shb - > ${id%_*}.bam
done


# 1.3 The results of Kallisto
# Generate three files：abundances.h5,abudances.tsv,run_info.json
# 1-abundance.h5
# HDF5 Binary format file, including run log information, expression abundance estimate, bootstrap estimate and transcript length information. The file can be directly read and processed by sleuth, or it can be converted into a plain text tsv format file using the kallisto h5dump command.
# 2-abundance.tsv
# Contains the pure text tsv format file with the header, the header is: target_id, length, eff_length, est_counts, tpm
# 3-run_info.json
# A log file in json format



### 2 Differential expression analysis
# 2.1 Load package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximport") 
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
library(tximport)
library(DEseq2)
library(limma)
library(edgeR)



# 2.2 Data format conversion. 
# extract chrA chrB chrD
grep 'TraesCS.A' Triticum_AC1-rep1/abundance.tsv > Triticum_AC1-rep1.chrA.counts.tsv
grep 'TraesCS.B' Triticum_AC1-rep1/abundance.tsv > Triticum_AC1-rep1.chrB.counts.tsv
grep 'TraesCS.D' Triticum_AC1-rep1/abundance.tsv > Triticum_AC1-rep1.chrD.counts.tsv
grep 'TraesCS.A' Triticum_AC1-rep2/abundance.tsv > Triticum_AC1-rep2.chrA.counts.tsv
grep 'TraesCS.B' Triticum_AC1-rep2/abundance.tsv > Triticum_AC1-rep2.chrB.counts.tsv
grep 'TraesCS.D' Triticum_AC1-rep2/abundance.tsv > Triticum_AC1-rep2.chrD.counts.tsv

# homeolog file
cut -f 1 homeolog.ABD.list > homeolog.A.list
awk ’{print $2"\t"$1;}’ homeolog.ABD.list > homeolog.B.list
awk ’{print $3"\t"$1;}’ homeolog.ABD.list > homeolog.D.list
# delete "transcript" in the file of homeolog:
sed s/transcript://g homeolog.A.list > homeolog.A.1.list
sed s/transcript://g homeolog.B.list > homeolog.B.1.list
sed s/transcript://g homeolog.D.list > homeolog.D.1.list

# homeolog-A-counts file
ls Triticum_AC1-rep*.chrA.counts.tsv |while read id
do
python ./tablize.py -a homeolog.A.1.list $id | sort -k1 > ${id%counts*}homeolog.tsv
done

# homeolog-B-counts file
ls Triticum_AC1-rep*.chrB.counts.tsv |while read id
do
python ./tablize.py -a homeolog.B.1.list $id | cut -f 2,3- | sort -k1 > ${id%counts*}homeolog.tsv
done

# homeolog-D-counts file
ls Triticum_AC1-rep*.chrD.counts.tsv |while read id
do
python ./tablize.py -a homeolog.D.1.list $id | cut -f 2,3- | sort -k1 > ${id%counts*}homeolog.tsv
done

# Check if the *homeolog.tsv format is correct, and then add the following column names
# target_id	length	eff_length	est_counts	tpm


# 2.3 Prepare differential expression data
# Copy the prepared *homeolog.tsv data to the counts folder
mkdir counts
cd counts
cp *homeolog.tsv ./

sample_id <- dir(file.path("/software/kallisto/DEG_ABD", "counts"))
sample_id 
kal_dirs <- file.path("/software/kallisto/DEG_ABD", "counts", sample_id)
kal_dirs

t2g <- read.table(file.path("/software/kallisto/test", "transcript_gene_Triticum_relation.txt"), header = TRUE)
t2g <- dplyr::rename(t2g,target_id = transcript_id, gene_id = gene_id)

library(tximport)
txi <- tximport(files = kal_dirs, type = "kallisto", tx2gene = t2g)
names(txi)
# [1] "abundance"           "counts"              "length"             
# [4] "countsFromAbundance"
head(txi$counts)                 # The order of the columns is based on the order of the files in the folder. If there are duplicates, the duplicate files are sorted together.



# 2.3.1 DEseq2
### (1) Grouping
sampleTable <- data.frame(condition = factor(rep(c("chrA", "chrB", "chrD"), each = 2)))
condition <- factor(c(rep("chrA",2), rep("chrB",2), rep("chrD",2)))
### (2) Matrix
rownames(sampleTable) <- colnames(txi$counts)
sampleTable                 # The row name of sampleTable corresponds to the column name of txi$counts
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# dds is now ready for DESeq() see DESeq2 vignette


# (3) Analysis of differential expression. 
# (3.1) difference analysis
if(T){
  suppressMessages(library(DESeq2)) 
  dds <- DESeq(dds)
  pdf("qc_dispersions.pdf")
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld) 
  write.csv(exprMatrix_rlog,'exprMatrix.rlog.csv' )
  
  normalizedCounts1 <- t( t(counts(dds)) / sizeFactors(dds) )
  # normalizedCounts2 <- counts(dds, normalized=T) # it's the same for the tpm value
  # we also can try cpm or rpkm from edgeR pacage
  exprMatrix_rpm=as.data.frame(normalizedCounts1) 
  head(exprMatrix_rpm)
  write.csv(exprMatrix_rpm,'exprMatrix.rpm.csv' )
  
  pdf("DEseq_RAWvsNORM.pdf")
  par(cex = 0.7)
  n.sample=ncol(txi$counts)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(txi$counts, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(txi$counts))
  hist(exprMatrix_rlog)
  dev.off()
  
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
  cor(as.matrix(txi$counts))
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(exprMatrix_rlog)))
  #install.packages("gplots",repos = "http://cran.us.r-project.org")
  library(gplots)
  pdf("qc-heatmap-samples.pdf")
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  cor(exprMatrix_rlog) 
  
  res <- results(dds, contrast=c("condition","chrA","chrB"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrA_chrB=as.data.frame(resOrdered)
  write.csv(DEG_chrA_chrB,"DEG_chrA_chrB_deseq2.results.csv")
  
  res <- results(dds, contrast=c("condition","chrA","chrD"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrA_chrD=as.data.frame(resOrdered)
  write.csv(DEG_chrA_chrD,"DEG_chrA_chrD_deseq2.results.csv")

  res <- results(dds, contrast=c("condition","chrB","chrD"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrB_chrD=as.data.frame(resOrdered)
  write.csv(DEG_chrB_chrD,"DEG_chrB_chrD_deseq2.results.csv")
}

### (3.2) Extract the difference matrix with significant genetic differences
# # chrA vs chrB
allDEG2 <- read.csv('DEG_chrA_chrB_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrA_chrB= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrA_chrB = diff_signif_DEG_chrA_chrB[order(diff_signif_DEG_chrA_chrB$log2FoldChange),]
diff_signif_DEG_chrA_chrB=na.omit(diff_signif_DEG_chrA_chrB)                       #Remove the line containing NA
write.csv(diff_signif_DEG_chrA_chrB,"DEG_chrA_chrB_deseq2_diff_signif.results.csv")

# chrA vs chrD
allDEG2 <- read.csv('DEG_chrA_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrA_chrD= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrA_chrD = diff_signif_DEG_chrA_chrD[order(diff_signif_DEG_chrA_chrD$log2FoldChange),]
diff_signif_DEG_chrA_chrD=na.omit(diff_signif_DEG_chrA_chrD)                       #Remove the line containing NA
write.csv(diff_signif_DEG_chrA_chrD,"DEG_chrA_chrD_deseq2_diff_signif.results.csv")

# chrB vs chrD
allDEG2 <- read.csv('DEG_chrB_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrB_chrD= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrB_chrD = diff_signif_DEG_chrB_chrD[order(diff_signif_DEG_chrB_chrD$log2FoldChange),]
diff_signif_DEG_chrB_chrD=na.omit(diff_signif_DEG_chrB_chrD)                       #Remove the line containing NA
write.csv(diff_signif_DEG_chrB_chrD,"DEG_chrB_chrD_deseq2_diff_signif.results.csv")

### (3.3) Comparison of the differentially expressed genes in screening results of chrA_chrB, chrA_chrD, chrB_chrD 
# Load chrA_chrB data
chrA_chrB <- read.csv('DEG_chrA_chrB_deseq2_diff_signif.results.csv',header = T , stringsAsFactors = F)
dim(chrA_chrB)                        # View matrix dimensions: several rows and columns # 62  8
A_vs_B = rownames(chrA_chrB)
A_vs_B = chrA_chrB[,2]               # Take A_vs_B as the second column, the gene name
# Load chrA_chrD data
chrA_chrD <- read.csv('DEG_chrA_chrD_deseq2_diff_signif.results.csv',header = T , stringsAsFactors = F)
dim(chrA_chrD)                        # 51   8
A_vs_D = rownames(chrA_chrD)
A_vs_D = chrA_chrD[,2]               
# Load chrB_chrD data
chrB_chrD <- read.csv('DEG_chrB_chrD_deseq2_diff_signif.results.csv',header = T , stringsAsFactors = F)
dim(chrB_chrD)                         # 60  8
B_vs_D = rownames(chrB_chrD)
B_vs_D = chrB_chrD[,2]              

venn.diagram(
  x = list(
    'A_vs_B(62)' = A_vs_B,
    'A_vs_D(51)' = A_vs_D,
    'B_vs_D(60)' = B_vs_D
  ),
  filename = 'VN.pdf',
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main = "chrA_chrB_chrD",
  main.cex = 1.2
)


### (3.4) Volcano Map
# chrA vs chrB
## Load the result file generated in DESeq2
# rm(list=ls())
resdata <- read.csv('DEG_chrA_chrB_deseq2.results.csv',header = T , stringsAsFactors = F)
threshold <- as.factor(ifelse(resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 2 ,
                              ifelse(resdata$log2FoldChange >= 2 ,'Up','Down'),'Not'))
pdf('DEG_chrA_chrB_deseq2.pdf')
ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),colour=threshold)) +
  xlab("log2(Fold Change)")+ylab("-log10(qvalue)") +
  geom_point(size = 0.5,alpha=1) +
  ylim(0,20) + xlim(-12,12) +
  scale_color_manual(values=c("green","grey", "red"))
dev.off()

# chrA vs chrD
resdata <- read.csv('DEG_chrA_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
threshold <- as.factor(ifelse(resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 2 ,
                              ifelse(resdata$log2FoldChange >= 2 ,'Up','Down'),'Not'))
pdf('DEG_chrA_chrD_deseq2.pdf')
ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),colour=threshold)) +
  xlab("log2(Fold Change)")+ylab("-log10(qvalue)") +
  geom_point(size = 0.5,alpha=1) +
  ylim(0,20) + xlim(-12,12) +
  scale_color_manual(values=c("green","grey", "red"))
dev.off()

# chrB vs chrD
resdata <- read.csv('DEG_chrB_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
threshold <- as.factor(ifelse(resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 2 ,
                              ifelse(resdata$log2FoldChange >= 2 ,'Up','Down'),'Not'))
pdf('DEG_chrB_chrD_deseq2.pdf')
ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),colour=threshold)) +
  xlab("log2(Fold Change)")+ylab("-log10(qvalue)") +
  geom_point(size = 0.5,alpha=1) +
  ylim(0,20) + xlim(-12,12) +
  scale_color_manual(values=c("green","grey", "red"))
dev.off()



# 2.3.2 edgeR
### (1) Grouping
sampleTable <- data.frame(condition = factor(rep(c("chrA", "chrB", "chrD"), each = 2)))
condition <- factor(c(rep("chrA",2), rep("chrB",2), rep("chrD",2)))
### (2) Matrix
cts <- txi$counts
normMat <- txi$length
# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide
### (3) Differential expression analysis
if(T){
  library(edgeR)
  y <- calcNormFactors(y)
  y$samples
  
  ## The calcNormFactors function normalizes for RNA composition by finding a set of scaling--TMM算法标准化
  ## factors for the library sizes that minimize the log-fold changes between the samples for most
  ## genes. The default method for computing these scale factors uses a trimmed mean of Mvalues
  ## (TMM) between each pair of samples
  
  pdf('edgeR_MDS.pdf')
  plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
  dev.off()
  
  # The glm approach to multiple groups is similar to the classic approach, but permits more general comparisons to be made
  
  dge=y
 
  design <- model.matrix(~condition, data = condition)
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(condition))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  
  lrt <- glmLRT(fit,  contrast=c(1,-1,0))     # chrA vs chrB
  nrDEG=topTags(lrt, n=nrow(txi$counts))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrA_chrB_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrA_chrB_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()

  lrt <- glmLRT(fit, contrast=c(1,0,-1) )      # chrA vs chrD
  nrDEG=topTags(lrt, n=nrow(txi$counts))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrA_chrD_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrA_chrD_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  lrt <- glmLRT(fit, contrast=c(0,1,-1) )      # chrB vs chrD
  nrDEG=topTags(lrt, n=nrow(txi$counts))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrB_chrD_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrB_chrD_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()
}



# 2.3.3 limma_voom
### (1) Grouping
sampleTable <- data.frame(condition = factor(rep(c("chrA", "chrB", "chrD"), each = 2)))
condition <- factor(rep(c("chrA", "chrB", "chrD"), each = 2))
### (2) Matrix
txi <- tximport(files = kal_dirs, type = "kallisto", tx2gene = t2g, countsFromAbundance = "lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)
# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
y <- calcNormFactors(y)
design <- model.matrix(~condition, data = condition)
v <- voom(y, design)
# v is now ready for lmFit() see limma User's Guide


### (3) Differential expression analysis
suppressMessages(library(limma))
colnames(design)=levels(factor(condition))
rownames(design)=colnames(txi$counts)

fit <- lmFit(v, design)

condition
cont.matrix=makeContrasts(contrasts=c('chrA-chrB','chrA-chrD','chrB-chrD'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='chrA-chrB', n=Inf)
DEG_chrA.vs.chrB_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrA.vs.chrB_limma_voom,"DEG_chrA-chrB_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='chrA-chrD', n=Inf)
DEG_chrA.vs.chrD_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrA.vs.chrD_limma_voom,"DEG_chrA-chrD_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='chrB-chrD', n=Inf)
DEG_chrB.vs.chrD_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrB.vs.chrD_limma_voom,"DEG_chrB-chrD_limma_voom.csv",quote = F)

pdf("limma_voom_RAWvsNORM.pdf")
txi$counts_new=v$E
par(cex = 0.7)
n.sample=ncol(txi$counts)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(txi$counts, col = cols,main="expression value",las=2)
boxplot(txi$counts_new, col = cols,main="expression value",las=2)
hist(as.matrix(txi$counts))
hist(txi$counts_new)
dev.off()


