
### Install software
# gffread
# https://github.com/gpertea/gffread.git
git clone https://github.com/gpertea/gffread.git
cd gffread
make release
echo "PATH=/home/minmin_li/software/gffread:$PATH" >> ~/.bashrc
source ~/.bashrc
gffread -h

# last
conda activate rnaseq
conda install -c bioconda last=1060 -y

# subread
conda activate rnaseq
conda install -c bioconda subread=2.0.0

# samtools
conda activate rnaseq
conda install -c bioconda samtools=1.3.1



### 1 Homeolog identification
CPU=8
REF=Triticum_aestivum.IWGSC.dna.toplevel
GTF=Triticum_aestivum.IWGSC.46
# Extract transcript sequences

# gff to gtf
gffread -T -o refseq.gtf $GTF.gff3
# Generate transcript file
gffread -g $REF.fa -w refseq.fa $GTF.gff3

# Generate refseq.chrA/B/D.gtf
grep '^.A' refseq.gtf > refseq.chrA.gtf
grep '^.B' refseq.gtf > refseq.chrB.gtf
grep '^.D' refseq.gtf > refseq.chrD.gtf
# or
# cat refseq.gtf | awk '($1=="1A" || $1=="2A" || $1=="3A" || $1=="4A" || $1=="5A" || $1=="6A" || $1=="7A"){ print}' > refseq.chrA.gtf
# cat refseq.gtf | awk '($1=="1B" || $1=="2B" || $1=="3B" || $1=="4B" || $1=="5B" || $1=="6B" || $1=="7B"){ print}' > refseq.chrB.gtf
# cat refseq.gtf | awk '($1=="1D" || $1=="2D" || $1=="3D" || $1=="4D" || $1=="5D" || $1=="6D" || $1=="7D"){ print}' > refseq.chrD.gtf

# Generate chrA/B/D.cds.fa
gffread -g $REF.fa -x chrA.cds.fa refseq.chrA.gtf
gffread -g $REF.fa -x chrB.cds.fa refseq.chrB.gtf
gffread -g $REF.fa -x chrD.cds.fa refseq.chrD.gtf

lastdb -uNEAR -R01 chrA_db chrA.cds.fa
lastdb -uNEAR -R01 chrB_db chrB.cds.fa
lastdb -uNEAR -R01 chrD_db chrD.cds.fa

lastal chrA_db -P$CPU -D10000000000 chrB.cds.fa | last-map-probs -m 0.49 > A.B.maf
lastal chrB_db -P$CPU -D10000000000 chrA.cds.fa | last-map-probs -m 0.49 > B.A.maf

lastal chrB_db -P$CPU -D10000000000 chrD.cds.fa | last-map-probs -m 0.49 > B.D.maf
lastal chrD_db -P$CPU -D10000000000 chrB.cds.fa | last-map-probs -m 0.49 > D.B.maf

lastal chrA_db -P$CPU -D10000000000 chrD.cds.fa | last-map-probs -m 0.49 > A.D.maf
lastal chrD_db -P$CPU -D10000000000 chrA.cds.fa | last-map-probs -m 0.49 > D.A.maf

python scripts/homeolog_genotypes.py -o A.vs.B -f exon -g refseq.gtf A.B.maf B.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o B.vs.A -f exon -g refseq.gtf B.A.maf A.B.maf # coordinates based on B

python scripts/homeolog_genotypes.py -o B.vs.D -f exon -g refseq.gtf B.D.maf D.B.maf # coordinates based on B
python scripts/homeolog_genotypes.py -o D.vs.B -f exon -g refseq.gtf D.B.maf B.D.maf # coordinates based on D

python scripts/homeolog_genotypes.py -o A.vs.D -f exon -g refseq.gtf A.D.maf D.A.maf # coordinates based on A
python scripts/homeolog_genotypes.py -o D.vs.A -f exon -g refseq.gtf D.A.maf A.D.maf # coordinates based on D

# Triple copy homeologs
perl triple_homeolog.pl A.vs.B.reciprocal_best B.vs.D.reciprocal_best A.vs.D.reciprocal_best > homeolog.ABD.list
# Subgenome unique transcripts
cat A.vs.B.reciprocal_best A.vs.D.reciprocal_best | cut -f1 | sort | uniq > A.vs.all.list
cat B.vs.A.reciprocal_best B.vs.D.reciprocal_best | cut -f1 | sort | uniq > B.vs.all.list
cat D.vs.A.reciprocal_best D.vs.B.reciprocal_best | cut -f1 | sort | uniq > D.vs.all.list
python scripts/tablize.py -v0 A.vs.all.list chrA.cds.fa > chrA.only.list
python scripts/tablize.py -v0 B.vs.all.list chrB.cds.fa > chrB.only.list
python scripts/tablize.py -v0 D.vs.all.list chrD.cds.fa > chrD.only.list



### 2 last
# Options: Genomic files are too large, with 1A 1B 1D genomic data.
cat Triticum_aestivum.IWGSC.dna.chr1A.fa Triticum_aestivum.IWGSC.dna.chr1A.fa Triticum_aestivum.IWGSC.dna.chr1A.fa > Triticum_aestivum.IWGSC.dna.chr1A1BD1.fa
lastdb8 -W 3 -uNEAR -R01 Taes Triticum_aestivum.IWGSC.dna.chr1A1BD1.fa
samtools dict -o Triticum_aestivum.IWGSC.dna.chr1A1BD1.dict Triticum_aestivum.IWGSC.dna.chr1A1BD1.fa

# Unzipt the sequencing fastq.gz file.
gzip -dc Triticum_AC1-rep1_1.fastq.gz | paste - - - - | perl -ne 'chomp; @t=split(/\s+/); @u=split(/\t/); print "$t[0]_1\n$u[1]\n+\n$u[3]\n";' > Triticum_AC1-rep1_1.fastq
gzip -dc Triticum_AC1-rep1_2.fastq.gz | paste - - - - | perl -ne 'chomp; @t=split(/\s+/); @u=split(/\t/); print "$t[0]_1\n$u[1]\n+\n$u[3]\n";' > Triticum_AC1-rep1_2.fastq

lastal8 -Q1 -D100 -P8 Taes Triticum_AC1-rep1_1.fastq Triticum_AC1-rep1_2.fastq | last-split8 -c 0 -t 0.004 -d 2 -m 1 -g Taes > ./Triticum_AC1-rep1.maf

maf-convert -f Triticum_aestivum.IWGSC.dna.chr1A.dict sam ./Triticum_AC1-rep1.maf | samtools view -Shb - > ./Triticum_AC1-rep1.maf.bam
samtools sort -o ./Triticum_AC1-rep1.maf.refsort.bam ./Triticum_AC1-rep1.maf.bam
samtools index -c ./Triticum_AC1-rep1.maf.refsort.bam

featureCounts -T 8 -Q 20 -t exon -g transcript_id -a refseq.chrA.gtf -o ./Triticum_AC1-rep1.chrA.counts.txt ./Triticum_AC1-rep1.maf.refsort.bam
featureCounts -T 8 -Q 20 -t exon -g transcript_id -a refseq.chrB.gtf -o ./Triticum_AC1-rep1.chrB.counts.txt ./Triticum_AC1-rep1.maf.refsort.bam
featureCounts -T 8 -Q 20 -t exon -g transcript_id -a refseq.chrD.gtf -o ./Triticum_AC1-rep1.chrD.counts.txt ./Triticum_AC1-rep1.maf.refsort.bam

python ../tablize.py -skip 1 -a -i 0 -c 6 ./*.chrA.counts.txt > Triticum_AC1-rep1.chrA.tsv
python ../tablize.py -skip 1 -a -i 0 -c 6 ./*.chrB.counts.txt > Triticum_AC1-rep1.chrB.tsv
python ../tablize.py -skip 1 -a -i 0 -c 6 ./*.chrD.counts.txt > Triticum_AC1-rep1.chrD.tsv


# In terms of chrA gene id
cut -f 1 homeolog.ABD.list > homeolog.A.list
python scripts/tablize.py -a homeolog.A.list Triticum_AC1-rep1.chrA.tsv | sort -k1 > last.chrA.homeolog.tsv

awk '{print $2"\t"$1;}' homeolog.ABD.list > homeolog.B.list
python scripts/tablize.py -a homeolog.B.list Triticum_AC1-rep1.chrB.tsv | cut -f 2,3- | sort -k1 > last.chrB.homeolog.tsv

awk '{print $3"\t"$1;}' homeolog.ABD.list > homeolog.D.list
python scripts/tablize.py -a homeolog.D.list Triticum_AC1-rep1.chrD.tsv | cut -f 2,3- | sort -k1 > last.chrD.homeolog.tsv



### 3 Analysis of the expression of differences.
# (1) Data format conversion.
sed "s/:/\t/g; s/transcript/chrA/g" last.chrA.homeolog.tsv > last.chrA.homeolog-1.tsv
sed "s/:/\t/g; s/transcript/chrB/g" last.chrB.homeolog.tsv > last.chrB.homeolog-1.tsv
sed "s/:/\t/g; s/transcript/chrD/g" last.chrD.homeolog.tsv > last.chrD.homeolog-1.tsv
cat last.chrA.homeolog-1.tsv last.chrB.homeolog-1.tsv last.chrD.homeolog-1.tsv > last.homeolog.txt

# (2) Load the package.
conda activate rnaseq
R
# View installed packages (not necessarily full)
(.packages()) 
# Update package (not necessarily updated)
update.packages("stats4")
# Install the package.
install.packages("reshape2")
install.packages("gplots")
install.packages("ggplot2")
install.packages("VennDiagram")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
# Install the dependency package.
BiocManager::install("DESeq2")
BiocManager::install("stats4")
BiocManager::install("BiocGenerics")
BiocManager::install("parallel")
BiocManager::install("BiocGenerics")
BiocManager::install("IRanges")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomeInfoDb")
BiocManager::install("SummarizedExperiment")
BiocManager::install("Biobase")
BiocManager::install("DelayedArray")
BiocManager::install("matrixStats")
BiocManager::install("BiocParallel")
BiocManager::install("DelayedArray")
# Load the dependency package.
rm(list=ls())
library(limma)
library(DESeq2)
library(stats4)
library(BiocGenerics)
library(parallel)
library(BiocGenerics)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(Biobase)
library(DelayedArray)
library(matrixStats)
library(BiocParallel)
library(DelayedArray)
library(grid)
library(futile.logger)

# Load the package used for the analysis.
library(reshape2)
library(edgeR)
library(DESeq2)
library(gplots)
library(ggplot2)
library(VennDiagram)

# (3) Prepare exprSet matrix for difference analysis
# DEseq2 edgeR limma analysis requires duplicate data
setwd("software/RNA_seq/20200518-rnaseq/LAST/1_result/DEseq")
a=read.table('last.homeolog-DEGs.txt',stringsAsFactors = F)
######################################################################
#chrA-1	TraesCS1A02G001800.1	2
#chrA-1	TraesCS1A02G002000.1	4
#chrA-1	TraesCS1A02G002200.1	25
#chrA-1	TraesCS1A02G002300.1	2
#chrA-1	TraesCS1A02G002500.1	14
#chrA-1	TraesCS1A02G003900.1	0
#chrA-2	TraesCS1A02G001800.1	20
#chrA-2	TraesCS1A02G002000.1	1
#chrA-2	TraesCS1A02G002200.1	29
#chrA-2	TraesCS1A02G002300.1	2
#chrA-2	TraesCS1A02G002500.1	14
#chrA-2	TraesCS1A02G003900.1	0
#chrB-1	TraesCS1A02G001800.1	2
#chrB-1	TraesCS1A02G002000.1	0
#chrB-1	TraesCS1A02G002200.1	40
#chrB-1	TraesCS1A02G002300.1	0
#chrB-1	TraesCS1A02G002500.1	1
#chrB-1	TraesCS1A02G003900.1	0
#chrB-2	TraesCS1A02G001800.1	2
#chrB-2	TraesCS1A02G002000.1	0
#chrB-2	TraesCS1A02G002200.1	40
#chrB-2	TraesCS1A02G002300.1	0
#chrB-2	TraesCS1A02G002500.1	1
#chrB-2	TraesCS1A02G003900.1	0
#chrD-1	TraesCS1A02G001800.1	0
#chrD-1	TraesCS1A02G002000.1	16
#chrD-1	TraesCS1A02G002200.1	0
#chrD-1	TraesCS1A02G002300.1	2
#chrD-1	TraesCS1A02G002500.1	0
#chrD-1	TraesCS1A02G003900.1	0
#chrD-2	TraesCS1A02G001800.1	0
#chrD-2	TraesCS1A02G002000.1	16
#chrD-2	TraesCS1A02G002200.1	0
#chrD-2	TraesCS1A02G002300.1	2
#chrD-2	TraesCS1A02G002500.1	0
#chrD-2	TraesCS1A02G003900.1	0
######################################################################

colnames(a)=c('sample','gene','reads')    # Add column names to the columns, sample is equivalent to different chromosomes chr
exprSet=dcast(a,gene~sample)
head(exprSet)
# write.table(exprSet,'stats.txt',quote=F,sep='\t')
 
geneLists=exprSet[,1]         # First column
exprSet=exprSet[,-1]          # Remove the first column
head(exprSet)

rownames(exprSet)=geneLists
write.csv(exprSet,'raw_reads_DEseq_matrix.csv') 

keepGene=rowSums(cpm(exprSet)>0) >=2
table(keepGene);dim(exprSet)
dim(exprSet[keepGene,])
exprSet=exprSet[keepGene,]
rownames(exprSet)=geneLists[keepGene]
str(exprSet)
group_list <- factor(c(rep("chrA",2), rep("chrB",2), rep("chrD",2)))
write.csv(exprSet,'filter_reads_DEGs_matrix.csv')


### (3.1) DEseq2
######################################################################
###################      Firstly for DEseq2      #####################
######################################################################
# Script reference URL:https://github.com/jmzeng1314/my-R/blob/master/10-RNA-seq-3-groups/hisat2_mm10_htseq.R
# Other related websites:
# http://www.360doc.com/content/18/0717/17/46931810_771191049.shtml
# https://www.jianshu.com/p/80c7bf8a2834
# https://www.jianshu.com/p/0dcc6030343e

### 1 Difference analysis
if(T){
  suppressMessages(library(DESeq2)) 
  (colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
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
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(group_list))])
  cor(as.matrix(exprSet))
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(exprMatrix_rlog)))
  #install.packages("gplots",repos = "http://cran.us.r-project.org")
  library(gplots)
  pdf("qc-heatmap-samples.pdf")
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[group_list], RowSideColors=mycols[group_list],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  cor(exprMatrix_rlog) 
  
  res <- results(dds, contrast=c("group_list","chrA","chrB"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrA_chrB=as.data.frame(resOrdered)
  write.csv(DEG_chrA_chrB,"DEG_chrA_chrB_deseq2.results.csv")
  
  res <- results(dds, contrast=c("group_list","chrA","chrD"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrA_chrD=as.data.frame(resOrdered)
  write.csv(DEG_chrA_chrD,"DEG_chrA_chrD_deseq2.results.csv")

  res <- results(dds, contrast=c("group_list","chrB","chrD"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_chrB_chrD=as.data.frame(resOrdered)
  write.csv(DEG_chrB_chrD,"DEG_chrB_chrD_deseq2.results.csv")
}


### 2 Extract the difference matrix with significant genetic differences
# # chrA vs chrB
allDEG2 <- read.csv('DEG_chrA_chrB_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrA_chrB= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrA_chrB = diff_signif_DEG_chrA_chrB[order(diff_signif_DEG_chrA_chrB$log2FoldChange),]
diff_signif_DEG_chrA_chrB=na.omit(diff_signif_DEG_chrA_chrB)                       # remove lines containing NA
write.csv(diff_signif_DEG_chrA_chrB,"DEG_chrA_chrB_deseq2_diff_signif.results.csv")

# chrA vs chrD
allDEG2 <- read.csv('DEG_chrA_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrA_chrD= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrA_chrD = diff_signif_DEG_chrA_chrD[order(diff_signif_DEG_chrA_chrD$log2FoldChange),]
diff_signif_DEG_chrA_chrD=na.omit(diff_signif_DEG_chrA_chrD)                       # remove lines containing NA
write.csv(diff_signif_DEG_chrA_chrD,"DEG_chrA_chrD_deseq2_diff_signif.results.csv")

# chrB vs chrD
allDEG2 <- read.csv('DEG_chrB_chrD_deseq2.results.csv',header = T , stringsAsFactors = F)
padj = 0.001
foldChange= 2
diff_signif_DEG_chrB_chrD= allDEG2[(allDEG2$padj < padj & 
                          (allDEG2$log2FoldChange>foldChange | allDEG2$log2FoldChange<(-foldChange))),]
diff_signif_DEG_chrB_chrD = diff_signif_DEG_chrB_chrD[order(diff_signif_DEG_chrB_chrD$log2FoldChange),]
diff_signif_DEG_chrB_chrD=na.omit(diff_signif_DEG_chrB_chrD)                       # remove lines containing NA
write.csv(diff_signif_DEG_chrB_chrD,"DEG_chrB_chrD_deseq2_diff_signif.results.csv")


### 3 Comparison of screening results of differentially expressed genes of chrA_chrB, chrA_chrD and chrB_chrD
# Loading chrA_chrB data
chrA_chrB <- read.csv('DEG_chrA_chrB_deseq2_diff_signif.results.csv',header = T , stringsAsFactors = F)
dim(chrA_chrB)                        # View matrix dimensions: several rows and columns # 62  8
A_vs_B = rownames(chrA_chrB)
A_vs_B = chrA_chrB[,2]               # take A_vs_B as the second column, that is, the gene name
# Loading chrA_chrD data
chrA_chrD <- read.csv('DEG_chrA_chrD_deseq2_diff_signif.results.csv',header = T , stringsAsFactors = F)
dim(chrA_chrD)                        # 51   8
A_vs_D = rownames(chrA_chrD)
A_vs_D = chrA_chrD[,2]               
# Loading chrB_chrD data
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


# 4 Volcano map
# chrA vs chrB
## loading the result file generated in DESeq2

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


### (3.2) edgeR
######################################################################
###################      Then  for edgeR        #####################
######################################################################
# https://github.com/jmzeng1314/my-R/blob/master/10-RNA-seq-3-groups/hisat2_mm10_htseq.R
# https://www.jianshu.com/p/5f94ae79f298?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation


if(T){
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  ## The calcNormFactors function normalizes for RNA composition by finding a set of scaling--TMM algorithm standardization.
  ## factors for the library sizes that minimize the log-fold changes between the samples for most
  ## genes. The default method for computing these scale factors uses a trimmed mean of Mvalues
  ## (TMM) between each pair of samples
  
  pdf('edgeR_MDS.pdf')
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  dev.off()
  
  # The glm approach to multiple groups is similar to the classic approach, but permits more general comparisons to be made
  
  dge=d
 
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  
  lrt <- glmLRT(fit,  contrast=c(1,-1,0))     # chrA vs chrB
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrA_chrB_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrA_chrB_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()

  lrt <- glmLRT(fit, contrast=c(1,0,-1) )      # chrA vs chrD
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrA_chrD_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrA_chrD_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  lrt <- glmLRT(fit, contrast=c(0,1,-1) )      # chrB vs chrD
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_chrB_chrD_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  pdf("DEG_chrB_chrD_edgeR.pdf")
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
  dev.off()
}




### (3.3) limma/voom
######################################################################
###################      Then  for limma/voom        #################
######################################################################
# https://github.com/jmzeng1314/my-R/blob/master/10-RNA-seq-3-groups/hisat2_mm10_htseq.R

suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('chrA-chrB','chrA-chrD','chrB-chrD'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
 
tempOutput = topTable(fit2, coef='chrA-chrB', n=Inf)
DEG_chrA-chrB_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrA-chrB_limma_voom,"DEG_chrA-chrB_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='chrA-chrD', n=Inf)
DEG_chrA-chrD_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrA-chrD_limma_voom,"DEG_chrA-chrD_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='chrB-chrD', n=Inf)
DEG_chrB-chrD_limma_voom = na.omit(tempOutput)
write.csv(DEG_chrB-chrD_limma_voom,"DEG_chrB-chrD_limma_voom.csv",quote = F)

pdf("limma_voom_RAWvsNORM.pdf")
exprSet_new=v$E
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprSet_new)
dev.off()



### Reference article：https://academic.oup.com/bib/article/21/2/395/5251019
### Reference script：https://github.com/tony-kuo/eagle/tree/master/scripts

