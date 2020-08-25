
### 1 Software installation
### sratoolkit
curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz
tar xvf sratoolkit.2.8.2-ubuntu64.tar.gz
cd sratoolkit.2.8.2-ubuntu64/bin      
echo "export PATH=$PATH:/software/sratoolkit.2.8.2-ubuntu64/bin" >> ~/.bashrc   #Add environment variables
source ~/.bashrc
### BWA
# https://sourceforge.net/projects/bio-bwa/files/
tar jxvf bwa-0.5.9rc1.tar.bz2
cd bwa-0.7.17/
make    # Software compilation
echo "export PATH=/software/bwa-0.7.17:$PATH" >> ~/.bashrc   #Add environment variables
source ~/.bashrc
### GSNAP-GMAP
# http://research-pub.gene.com/gmap/src/README
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz
tar xf gmap-gsnap-2018-07-04.tar.gz
cd gmap-2018-07-04/
./configure --prefix=$HOME/opt/biosoft/gmap
make -j 20
make install
### BamBam (PolyCat)
# https://sourceforge.net/p/bambam/activity/
# Carefully read the install and README in the installation package
# First install SAMtools (including HTSlib), BAMtools, and Zlib
# Add LIBRARY_PATH
echo "export LIBRARY_PATH=/home/minmin_li/programs/samtools/samtools-1.9:/software/RSEM-1.3.3/samtools-1.3/htslib-1.3:/software/eagle/htslib:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/home/minmin_li/miniconda3/envs/rnaseq/lib:/software/zlib-1.2.11:/home/minmin_li/programs/samtools/samtools-1.9/win32:/home/minmin_li/programs/velvet/v1.2.10/x86_64/third-party/zlib-1.2.3:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/home/minmin_li/miniconda3/envs/rnaseq/lib:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/software/bamtools/build/src/api:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/software/bamtools/build/src/utils:$LIBRARY_PATH" >> ~/.bashrc
source ~/.bashrc



### 2 Prepare the data to be analyzed
### download
#  HANDS: a tool for genome-wide discovery of subgenome- specific base-identity in polyploids. BMC Genomics 14: 653. doi:10.1186/1471-2164-14-653.
# [Download](https://www.ncbi.nlm.nih.gov/sra/?term=SRA097144&utm_source=gquery&utm_medium=search)
# SRR949820.sra
# SRR949821.sra
# SRR949822.sra
# SRR949823.sra
# Unzip
# --split-3: If it is paired-end sequencing data, output two files, if not, output only one file
# --gzip: The output format is gzip compressed file (fastqc software can directly recognize gzip compressed files)
# -O|--outdir: Output to specified folder
# -A: File name after decompression
fastq-dump /HANDS/SRR949820/SRR949820.sra --gzip --split-3 -O ../raw_data
### Extract part of the data for analysis
head -n 400000 SRR949820_1.fastq.gz > SRR949820.1_1.fastq.gz
head -n 400000 SRR949820_2.fastq.gz > SRR949820.1_2.fastq.gz
head -n 400000 SRR949821_1.fastq.gz > SRR949821.1_1.fastq.gz
head -n 400000 SRR949821_2.fastq.gz > SRR949821.1_2.fastq.gz
head -n 400000 SRR949822_1.fastq.gz > SRR949822.1_1.fastq.gz
head -n 400000 SRR949822_2.fastq.gz > SRR949822.1_2.fastq.gz
head -n 400000 SRR949823_1.fastq.gz > SRR949823.1_1.fastq.gz
head -n 400000 SRR949823_2.fastq.gz > SRR949823.1_2.fastq.gz



### 3 Data quality control filtering
### (1)fastqc
# Options:
# -o --outdir Output directory, you need to create the directory yourself
# -noextract Do not unzip files
# -f Specify the type of input file, support fastq|bam|sam files in three formats, and automatically recognize it by default
ls *.fastq.gz | while read id;
do
fastqc -o /1_clean_reads/fastqc_1/ -t 3 -noextract -f fastq $id
done
### (2) Excision joint sequence
# Options:
# PE: Filter pair-end sequencing data (if the data is single-end sequencing, use SE).
# phred33: The quality value format of the Fastq file is phred33. Generally, the format of the second-generation sequencing data is basically phred33. If you are not sure, you can consult the sequencing company.
# threads: Set the number of threads.
# seq*.fq.gz: Fastq files to be filtered.
# seq*.clean.fq.gz: Fastq file after filtering.
# ILLUMINACLIP: ./adapters/TruSeq3-PE.fa refers to the removal of the TruSeq3 adapter sequence under the Illumina sequencing platform. For specific adapter sequences, please consult the sequencing company. The 3 numbers (2:30:10) following the adaptor sequence indicate that when the adaptor sequence is aligned, two base mismatches are allowed. If the two reads of paired-end sequencing match the adaptor sequence more than 30% , It will be cut off, and if the match rate of a single read with the linker sequence exceeds 10%, it will also be cut off.
# SLIDINGWINDOW:4:15: It means that the 4bp window is used for sliding window statistics, and the sequence of the window whose average quality of excision bases is lower than 15 and beyond.
# LEADING:3: Remove leading low quality or N bases (lower quality 3).
# TRAILING:3: Remove the tailing low quality or N bases (lower than quality 3).
# MINLEN:36: Indicates that the reads whose length is less than 36 after filtering are removed.

# SRR949820   Triticum aestivum AABBDD 100bp
java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 SRR949820.1_1.fastq.gz SRR949820.1_2.fastq.gz ./SRR949820.1_1.clean.fastq.gz ./SRR949820.1_1.unpaired.fastq.gz ./SRR949820.1_2.clean.fastq.gz ./SRR949820.1_2.unpaired.fastq.gz ILLUMINACLIP:/software/trimmomatic/v0.39/x86_64/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:75

# SRR949821   Triticum urartu AA 51bp
java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 SRR949821.1_1.fastq.gz SRR949821.1_2.fastq.gz ./SRR949821.1_1.clean.fastq.gz ./SRR949821.1_1.unpaired.fastq.gz ./SRR949821.1_2.clean.fastq.gz ./SRR949821.1_2.unpaired.fastq.gz ILLUMINACLIP:/software/trimmomatic/v0.39/x86_64/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

# SRR949822   Aegilops speltoides BB 51bp
java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 SRR949822_1.fastq.gz SRR949822_2.fastq.gz SRR949822_1.clean.fq.gz SRR949822_1.unpaired.fq.gz SRR949822_2.clean.fq.gz SRR949822_2.unpaired.fq.gz ILLUMINACLIP:/software/trimmomatic/v0.39/x86_64/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

# SRR949823   Aegilops tauschii DD 51bp
java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 SRR949823_1.fastq.gz SRR949823_2.fastq.gz SRR949823_1.clean.fq.gz SRR949823_1.unpaired.fq.gz SRR949823_2.clean.fq.gz SRR949823_2.unpaired.fq.gz ILLUMINACLIP:/software/trimmomatic/v0.39/x86_64/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

### (3)fastqc
ls *.clean.fastq.gz | while read id;
do
fastqc -o /home/minmin_li/Bioinformatics/RNA_seq/20200518-rnaseq/HANDS/20200525/1_clean_reads/fastqc_2 -t 3 -noextract -f fastq $id
done



### 4 Build SNP_index
### (1) Alignment of diploid reads to genomes of Aegilops tauschii(DD)
# BWA build index 
# After the index is established, five files will be generated with the suffix of bwt pac ann amb sa
bwa index Aegilops_tauschii.Aet_v4.0.cdna.all.fa.gz
# BWA alignment
bwa mem Aegilops_tauschii.Aet_v4.0.cdna.all.fa.gz SRR949821.1_1.clean.fastq.gz SRR949821.1_2.clean.fastq.gz > SRR949821.1.sam
bwa mem Aegilops_tauschii.Aet_v4.0.cdna.all.fa.gz SRR949822.1_1.clean.fastq.gz SRR949822.1_2.clean.fastq.gz > SRR949822.1.sam
bwa mem Aegilops_tauschii.Aet_v4.0.cdna.all.fa.gz SRR949823.1_1.clean.fastq.gz SRR949823.1_2.clean.fastq.gz > SRR949823.1.sam
# sam--bam (Version: 0.1.20)
ls *.sam |while read id
do
samtools view -bS $id > ${id%.*}.bam
done
# bam--sort
samtools sort SRR949821.1.bam SRR949821.sorted
samtools sort SRR949822.1.bam SRR949822.sorted
samtools sort SRR949823.1.bam SRR949823.sorted
# bam--index
samtools index SRR949821.sorted.bam
samtools index SRR949822.sorted.bam
samtools index SRR949823.sorted.bam

### (2) GSNAP builds an index for genome files
mkdir genome
# -d indicate the name of the database
# gmap_build -d reference reference.fa -D ./reference
# Genome files do not support compressed files
/software/gsnap/v20200408/gmap-2020-04-08/util/gmap_build -d genome Aegilops_tauschii.Aet_v4.0.cdna.all.fa -D ./genome

### (2) Build SNP_index
/software/bambam-1.4/bambam/bin/interSnp SRR949821.sorted.bam SRR949822.sorted.bam SRR949823.sorted.bam -t 5 > snps.txt
/software/bambam-1.4/bambam/scripts/snp2gsnap.pl snps.txt Aegilops_tauschii.Aet_v4.0.cdna.all.fa | /software/iit_store -o new_snps.iit
mkdir snp-index
/software/snpindex -D ../0_data/genome -d genome -V ./snp-index -v new_snps new_snps.iit



### 5 GSNAP SNP-tolerant alignment
# Alignment of hexaploid reads to diploid genomes
/software/gsnap --gunzip -n 1 -N 1 -Q -t 10 --merge-distant-samechr --use-sarray=0 -D ../0_data/genome -d genome -V ./snp-index -v new_snps -A sam SRR949820.1_1.clean.fastq.gz SRR949820.1_2.clean.fastq.gz > SRR949820.1.sam

# sam-bam-sort-index (1.3.1 version samtools-conda)
samtools view -Sb -@ 10 SRR949820.1.sam -o SRR949820.1.bam
samtools sort SRR949820.1.bam > SRR949820.1.sorted.bam
samtools index SRR949820.1.sorted.bam



### 6 PolyCat alignment
# -p 1 for pairend
/software/bambam-1.4/bambam/bin/polyCat -x 1 -p 1 -s ./snp-index SRR949820.1.sorted.bam


