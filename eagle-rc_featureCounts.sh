# Install software
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

# star
conda activate rnaseq
conda install -c bioconda star=2.7.3a -y

# eagle
git clone https://github.com/tony-kuo/eagle.git
cd eagle
git clone https://github.com/samtools/htslib.git
make



### 1 Homeolog identification
CPU=8
REF=Triticum_aestivum.IWGSC.dna.toplevel
GTF=Triticum_aestivum.IWGSC.46
# Extract transcript sequences

# gff to gtf
gffread -T -o refseq.gtf $GTF.gff3
# Create a genome index
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


### 2 Subgenome alignment
# optional: genome file and gtf file are too large, extract part of the data
# Genome data: extract the first 57555 lines of the 1A 1B 1D subgenome
head -n 57555 Triticum_aestivum.IWGSC.dna.chromosome.1A.fa > Triticum_aestivum.IWGSC.dna.chromosome.2A.fa
head -n 57555 Triticum_aestivum.IWGSC.dna.chromosome.1B.fa > Triticum_aestivum.IWGSC.dna.chromosome.2B.fa
head -n 57555 Triticum_aestivum.IWGSC.dna.chromosome.1D.fa > Triticum_aestivum.IWGSC.dna.chromosome.2D.fa
# gtf file: extract the first 1161 lines of 1A 1D subgenomic gtf file and 1B subgenomic gtf file
grep '^1A' refseq.gtf > refseq.chr1A.gtf
grep '^1B' refseq.gtf > refseq.chr1B.gtf
grep '^1D' refseq.gtf > refseq.chr1D.gtf
head -n 1161 refseq.chr1A.gtf > refseq.chr2A.gtf
head -n 1161 refseq.chr1B.gtf > refseq.chr2B.gtf
head -n 900 refseq.chr1D.gtf > refseq.chr2D.gtf

# Build index
# --genomeDir ：index Output path
mkdir Taes_chr2A Taes_chr2B Taes_chr2D
STAR --runMode genomeGenerate --genomeDir Taes_chr2A --genomeFastaFiles Triticum_aestivum.IWGSC.dna.chromosome.2A.fa --sjdbGTFfile refseq.chr2A.gtf --runThreadN 1
STAR --runMode genomeGenerate --genomeDir Taes_chr2B --genomeFastaFiles Triticum_aestivum.IWGSC.dna.chromosome.2B.fa --sjdbGTFfile refseq.chr2B.gtf --runThreadN 1
STAR --runMode genomeGenerate --genomeDir Taes_chr2D --genomeFastaFiles Triticum_aestivum.IWGSC.dna.chromosome.2D.fa --sjdbGTFfile refseq.chr2D.gtf --runThreadN 1

# alignment analysis
STAR --genomeDir Taes_chr2A --readFilesCommand zcat --readFilesIn Triticum_AC1-rep1_1.fastq.gz Triticum_AC1-rep1_2.fastq.gz \
--outFileNamePrefix star_Triticum_AC1-rep1-chr2A --runThreadN 1 --genomeLoad NoSharedMemory \
--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
--outStd SAM | samtools view -Shb - > Triticum_AC1-rep1-chr2A.bam

STAR --genomeDir Taes_chr2B --readFilesCommand zcat --readFilesIn Triticum_AC1-rep1_1.fastq.gz Triticum_AC1-rep1_2.fastq.gz \
--outFileNamePrefix star_Triticum_AC1-rep1-chr2B --runThreadN 1 --genomeLoad NoSharedMemory \
--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
--outStd SAM | samtools view -Shb - > Triticum_AC1-rep1-chr2B.bam

STAR --genomeDir Taes_chr2D --readFilesCommand zcat --readFilesIn Triticum_AC1-rep1_1.fastq.gz Triticum_AC1-rep1_2.fastq.gz \
--outFileNamePrefix star_Triticum_AC1-rep1-chr2D --runThreadN 1 --genomeLoad NoSharedMemory \
--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout \
--outStd SAM | samtools view -Shb - > Triticum_AC1-rep1-chr2D.bam

# Sort index for bam file
samtools sort -o ./Triticum_AC1-rep1-chr2A.refsort.bam ./Triticum_AC1-rep1-chr2A.bam
samtools index -c ./Triticum_AC1-rep1-chr2A.refsort.bam

samtools sort -o ./Triticum_AC1-rep1-chr2B.refsort.bam ./Triticum_AC1-rep1-chr2B.bam
samtools index -c ./Triticum_AC1-rep1-chr2B.refsort.bam

samtools sort -o ./Triticum_AC1-rep1-chr2D.refsort.bam ./Triticum_AC1-rep1-chr2D.bam
samtools index -c ./Triticum_AC1-rep1-chr2D.refsort.bam



### 3 EAGLE-RC
# Put the appropriate vcfs to the corresponding dir, i.e. A.vs.*.gtf.vcf in chrA
# Using Triticum_aestivum.IWGSC.dna.chromosome.2A.fa will prompt the lack of 3A 5A and other genomic data, re-use the whole genome data Triticum_aestivum.IWGSC.dna.toplevel.fa here
# eagle
eagle -t 8 -a Triticum_AC1-rep1-chr2A.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v A.vs.B.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2A.A.vs.B.txt 2> Triticum_AC1-rep1-chr2A.A.vs.B.readinfo.txt

eagle -t 8 -a Triticum_AC1-rep1-chr2B.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v B.vs.A.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2B.B.vs.A.txt 2> Triticum_AC1-rep1-chr2B.B.vs.A.readinfo.txt

eagle -t 8 -a Triticum_AC1-rep1-chr2A.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v A.vs.D.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2A.A.vs.D.txt 2> Triticum_AC1-rep1-chr2A.A.vs.D.readinfo.txt

eagle -t 8 -a Triticum_AC1-rep1-chr2D.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v D.vs.A.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2D.D.vs.A.txt 2> Triticum_AC1-rep1-chr2D.D.vs.A.readinfo.txt

eagle -t 8 -a Triticum_AC1-rep1-chr2D.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v D.vs.B.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2D.D.vs.B.txt 2> Triticum_AC1-rep1-chr2D.D.vs.B.readinfo.txt

eagle -t 8 -a Triticum_AC1-rep1-chr2B.refsort.bam -r /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/0_data/Triticum_aestivum.IWGSC.dna.toplevel.fa -v B.vs.D.gtf.vcf --splice --rc 1> Triticum_AC1-rep1-chr2B.B.vs.D.txt 2> Triticum_AC1-rep1-chr2B.B.vs.D.readinfo.txt

# eagle-rc
eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2A.A.vs.B -v Triticum_AC1-rep1-chr2A.A.vs.B.txt Triticum_AC1-rep1-chr2A.A.vs.B.readinfo.txt > Triticum_AC1-rep1-chr2A.A.vs.B.list

eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2A.A.vs.D -v Triticum_AC1-rep1-chr2A.A.vs.D.txt Triticum_AC1-rep1-chr2A.A.vs.D.readinfo.txt > Triticum_AC1-rep1-chr2A.A.vs.D.list

eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2B.B.vs.D -v Triticum_AC1-rep1-chr2B.B.vs.D.txt Triticum_AC1-rep1-chr2B.B.vs.D.readinfo.txt > Triticum_AC1-rep1-chr2B.B.vs.D.list

eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2B.B.vs.A -v Triticum_AC1-rep1-chr2B.B.vs.A.txt Triticum_AC1-rep1-chr2B.B.vs.A.readinfo.txt > Triticum_AC1-rep1-chr2B.B.vs.A.list

eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2D.D.vs.B -v Triticum_AC1-rep1-chr2D.D.vs.B.txt Triticum_AC1-rep1-chr2D.D.vs.B.readinfo.txt > Triticum_AC1-rep1-chr2D.D.vs.B.list

eagle-rc -a Triticum_AC1-rep1-chr2A.refsort.bam --listonly -o Triticum_AC1-rep1-chr2D.D.vs.A -v Triticum_AC1-rep1-chr2D.D.vs.A.txt Triticum_AC1-rep1-chr2D.D.vs.A.readinfo.txt > Triticum_AC1-rep1-chr2D.D.vs.A.list

# Alignment-quantification
mkdir eagle
# -o --out String  Prefix for output BAM files.
python /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/eagle-rc/ref3_consensus.py --pe -u -d -o eagle/Triticum_AC1-rep1.ref \
-A Triticum_AC1-rep1-chr2A.A.vs.B.list Triticum_AC1-rep1-chr2A.A.vs.D.list \
-B Triticum_AC1-rep1-chr2B.B.vs.A.list Triticum_AC1-rep1-chr2B.B.vs.D.list \
-D Triticum_AC1-rep1-chr2D.D.vs.A.list Triticum_AC1-rep1-chr2D.D.vs.B.list
# Generate Triticum_AC1-rep1.ref.chrA.list, Triticum_AC1-rep1.ref.chrB.list, Triticum_AC1-rep1.ref.chrD.list files

eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2A.refsort.bam -o eagle/Triticum_AC1-rep1.chrA eagle/Triticum_AC1-rep1.ref.chrA.list
eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2B.refsort.bam -o eagle/Triticum_AC1-rep1.chrB eagle/Triticum_AC1-rep1.ref.chrB.list
eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2D.refsort.bam -o eagle/Triticum_AC1-rep1.chrD eagle/Triticum_AC1-rep1.ref.chrD.list
# Generate Triticum_AC1-rep1.chrA.ref.bam, Triticum_AC1-rep1.chrB.ref.bam, Triticum_AC1-rep1.chrD.ref.bam files

featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2A.gtf -o eagle/Triticum_AC1-rep1-chr2A.counts.txt eagle/Triticum_AC1-rep1.chrA.ref.bam
featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2B.gtf -o eagle/Triticum_AC1-rep1-chr2B.counts.txt eagle/Triticum_AC1-rep1.chrB.ref.bam
featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2D.gtf -o eagle/Triticum_AC1-rep1-chr2D.counts.txt eagle/Triticum_AC1-rep1.chrD.ref.bam
# Generate Triticum_AC1-rep1-chr2A.counts.txt, Triticum_AC1-rep1-chr2B.counts.txt, Triticum_AC1-rep1-chr2D.counts.txt files

# Double and triple homeolog counts
cd eagle
python /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/eagle-rc/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2A.counts.txt > eagle.chr2A.tsv
python /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/eagle-rc/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2B.counts.txt > eagle.chr2B.tsv
python /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/eagle-rc/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2D.counts.txt > eagle.chr2D.tsv
# Generate eagle.chr2A.tsv, eagle.chr2B.tsv, and eagle.chr2D.tsv files

# Triple homeologs counts in terms of chrA gene id
# cd eagle
cut -f 1 homeolog.ABD.list > homeolog.A.list
python scripts/tablize.py -a homeolog.A.list eagle.chr2A.tsv | sort -k1 > eagle.chrA.homeolog.tsv
awk '{print $2"\t"$1;}' homeolog.ABD.list > homeolog.B.list
python scripts/tablize.py -a homeolog.B.list eagle.chr2B.tsv | cut -f 2,3- | sort -k1 > eagle.chrB.homeolog.tsv
awk '{print $3"\t"$1;}' homeolog.ABD.list > homeolog.D.list
python scripts/tablize.py -a homeolog.D.list eagle.chr2D.tsv | cut -f 2,3- | sort -k1 > eagle.chrD.homeolog.tsv
# Generated eagle.chrA.homeolog.tsv, eagle.chrB.homeolog.tsv, eagle.chrD.homeolog.tsv files, import DEseq2, edgeR, or limma_voom program package for gene differential expression analysis


# Subgenome unique mapped reads
# cd eagle
echo "" > dummy.txt
eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2A.refsort.bam -u Triticum_AC1-rep1-chr2B.refsort.bam,Triticum_AC1-rep1-chr2D.refsort.bam -o eagle/Triticum_AC1-rep1.chrA.only dummy.txt
eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2B.refsort.bam -u Triticum_AC1-rep1-chr2A.refsort.bam,Triticum_AC1-rep1-chr2D.refsort.bam -o eagle/Triticum_AC1-rep1.chrB.only dummy.txt
eagle-rc --refonly --readlist -a Triticum_AC1-rep1-chr2D.refsort.bam -u Triticum_AC1-rep1-chr2A.refsort.bam,Triticum_AC1-rep1-chr2B.refsort.bam -o eagle/Triticum_AC1-rep1.chrD.only dummy.txt


featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2A.gtf -o eagle/Triticum_AC1-rep1-chr2A.only.counts.txt eagle/Triticum_AC1-rep1.chrA.only.ref.bam
featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2B.gtf -o eagle/Triticum_AC1-rep1-chr2B.only.counts.txt eagle/Triticum_AC1-rep1.chrB.only.ref.bam
featureCounts -T 8 -t exon -g transcript_id -a refseq.chr2D.gtf -o eagle/Triticum_AC1-rep1-chr2D.only.counts.txt eagle/Triticum_AC1-rep1.chrD.only.ref.bam


python scripts/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2A.only.counts.txt > eagle.chrA.only.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2B.only.counts.txt > eagle.chrB.only.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 ./Triticum_AC1-rep1-chr2D.only.counts.txt > eagle.chrD.only.tsv


# Subgenome unique mapped reads in subgenome unique genes (i.e. non-homeologs)
# cd eagle
python scripts/tablize.py -a chrA.only.list eagle.chr2A.tsv > eagle.chrA.only.unique.tsv
python scripts/tablize.py -a chrB.only.list eagle.chr2B.tsv > eagle.chrB.only.unique.tsv
python scripts/tablize.py -a chrD.only.list eagle.chr2D.tsv > eagle.chrD.only.unique.tsv



### 4 Reference article：https://academic.oup.com/bib/article/21/2/395/5251019
### Script：https://github.com/tony-kuo/eagle/tree/master/scripts


























