### 1 Download and install
# http://deweylab.github.io/RSEM/
make
make install
# Compile EBSeq
make ebseq
echo "export PATH=/home/minmin_li/software/RSEM-1.3.3/:$PATH" >> ~/.bashrc
source ~/.bashrc



### 2 Build References
# RSEM can call bowtie, bowtie2, STAR three comparison tools; here bowtie2 is used
rsem-prepare-reference -gtf ./Triticum_aestivum.IWGSC.46.gtf --bowtie2 ./Triticum_aestivum.IWGSC.dna.fa ./index/Taes



### 3 Calculate expression
# Use RSEM's rsem-calculate-expression command to perform bowtie2 alignment of reads and quantify the expression level
# -append-names is used to append the gene name and transcript name to the result
# -output-genome-bam: output the bam file based on the gene level in the result (the default is only the bam file at the transcription level)

# rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
for j in $( ls fastq/ | grep '_1.fastq.gz' );	do
rsem-calculate-expression -p 8 --bowtie2 --time --append-names --output-genome-bam \
--paired-end fastq/$j fastq/${j%%_1.fastq.gz}_2.fastq.gz ./index/Taes ./${j%%_1.fastq.gz}
done



### 4 Output
# Each sample has two bam files, Triticum_AC1-rep1.genome.bam and Triticum_AC1-rep1.transcript.bam. This is because of the --output-genome-bam parameter; and you can also see that the two bam files are more common The bam files after the transcriptome comparison are much larger. This is because RSEM uses the -k 200 parameter when calling the bowtie2 software for comparison, so RSEM is different from some other software, which considers multiple comparison reads.
# Triticum_AC1-rep1.genes.results and Triticum_AC1-rep1.isoforms.results files are quantitative results based on gene level and transcription level, respectively



### 5 Data format conversion
# Extract chrA chrB chrD
# $(j%%.result*), remove the characters on the right of the rightmost .result
for j in $( ls ./ | grep 'isoforms.results' );	do
grep 'TraesCS.A' $j > ${j%%.results*}.chrA.results
done

for j in $( ls ./ | grep 'isoforms.results' );	do
grep 'TraesCS.B' $j > ${j%%.results*}.chrB.results
done

for j in $( ls ./ | grep 'isoforms.results' );	do
grep 'TraesCS.D' $j > ${j%%.results*}.chrD.results
done

for j in $( ls ./ | grep 'genes.results' );	do
grep 'TraesCS.A' $j > ${j%%.results*}.chrA.results
done

for j in $( ls ./ | grep 'genes.results' );	do
grep 'TraesCS.B' $j > ${j%%.results*}.chrB.results
done

for j in $( ls ./ | grep 'genes.results' );	do
grep 'TraesCS.D' $j > ${j%%.results*}.chrD.results
done

# homeolog file
cut -f 1 homeolog.ABD.list > homeolog.A.list
awk ’{print $2"\t"$1;}’ homeolog.ABD.list > homeolog.B.list
awk ’{print $3"\t"$1;}’ homeolog.ABD.list > homeolog.D.list
# delete transcript of the homeolog file:
sed s/transcript://g homeolog.A.list > homeolog.A.1.list
sed s/transcript://g homeolog.B.list > homeolog.B.1.list
sed s/transcript://g homeolog.D.list > homeolog.D.1.list

# homeolog-A-counts file
ls *isoforms.chrA.results |while read id
do
python ./tablize.py -a homeolog.A.1.list $id | sort -k1 > ${id%results*}homeolog.results
done

# homeolog-B-counts file
ls *isoforms.chrB.results |while read id
do
python ./tablize.py -a homeolog.B.1.list $id | cut -f 2,3- | sort -k1 > ${id%results*}homeolog.results
done

# homeolog-D-counts file
ls *isoforms.chrD.results |while read id
do
python ./tablize.py -a homeolog.D.1.list $id | cut -f 2,3- | sort -k1 > ${id%results*}homeolog.results
done



### 6 Differential Expression Analysis using EBSeq
# Reference URL: https://www.jianshu.com/p/3e96d0f8fe71
# analysis of isoforms level differences
# Options: If there are two consecutive tab keys in the file *isoforms.chrA.homeolog.results, remove them with the following command, otherwise the resulting matrix is incorrect.
cat Triticum_AC1-rep2.isoforms.chrA.homeolog.results | tr -s '\t' > Triticum_AC1-rep2.isoforms.chrA1.homeolog.results
cat Triticum_AC1-rep1.isoforms.chrA.homeolog.results | tr -s '\t' > Triticum_AC1-rep1.isoforms.chrA1.homeolog.results

# (1) Extract the columns in each sample isoforms quantitative results expected_count to form a data matrix.
rsem-generate-data-matrix \
Triticum_AC1-rep1.isoforms.chrA1.homeolog.results Triticum_AC1-rep2.isoforms.chrA1.homeolog.results \
Triticum_AC1-rep1.isoforms.chrB.homeolog.results Triticum_AC1-rep2.isoforms.chrB.homeolog.results \
Triticum_AC1-rep1.isoforms.chrD.homeolog.results Triticum_AC1-rep2.isoforms.chrD.homeolog.results \
> ./GeneMat.txt

# (2) run-ebseq call EBseq for inspection
# 2,2,2 means 3 conditions, each condition has two repetitions; the order should be the same as the order of the condition indicated by the input file in the previous step
rsem-run-ebseq GeneMat.txt 2,2,2 GeneMat.results 
# Will get 4 files: GeneMat.results.condmeans, GeneMat.results, GeneMat.results.pattern, GeneMat.results.normalized_data_matrix

### GeneMat.results.pattern  
# *******************************
# "C1"    "C2"    "C3"
# "Pattern1"      1       1       1
# "Pattern2"      1       1       2
# "Pattern3"      1       2       1
# "Pattern4"      1       2       2
# "Pattern5"      1       2       3
# *******************************
# In the case of Pattern4, 1 2 2 indicates a gene expression: C1 is different from C2, C2 is the same as C3, and if there are differences in gene expression in the three types of condition, these are the cases.

### GeneMat.results
# "Pattern1"      "Pattern2"      "Pattern3"      "Pattern4"      "Pattern5"      "MAP"   "PPDE"
# The first column is the name of each gene, and the next five columns are the probability that the gene will match the party.
# "MAP" is the most likely mode of the gene; "PPDE"：posterior probability of being differentially expressed，the bigger the better.


### GeneMat.results.condmeans
# For each sample to combine the quantitative results after repetition, this result can be used with other tools for differential testing.
# *************************************************************************************
# "C1"    "C2"    "C3"
# "TraesCS1A02G096500.2"  0       1.96561377486638        302.775750524496
# "TraesCS1A02G348000.1"  74.6957365143482        10.3138887875017        107.689533956856
# "TraesCS2A02G109900.1"  29.1225097688886        163.015999330274        50.5862871503378
# "TraesCS4A02G283500.1"  0       0       40.7985187674878
# "TraesCS7A02G316200.1"  2.47966359816966        68.4194392774585        63.6670707799793
# "TraesCS1A02G005700.1"  222.628223296114        17.6702202263544        0
# "TraesCS1A02G007500.1"  87.3293858601812        15.7228798241867        0
# "TraesCS1A02G089900.1"  0       41.2291602783308        0
# **************************************************************************************

### GeneMat.results.normalized_data_matrix
# The matrix after the averageization.

# (3) Control_fdr control FDR (error detection rate) to select differential genes.
# The records in the GeneMat.results file with PPDE greater than 0.95 are extracted.
# controlling the false discovery rate (FDR) at level 0.05--controls the error discovery rate to 0.05.
rsem-control-fdr GeneMat.results 0.05 GeneMat.de.txt
# According to GeneMat.results.pattern, extract two genes that are different.
cat GeneMat.de.txt | grep -E 'Pattern3|Pattern4|Pattern5' > DEG_chrA-chrB_EBSeq_rsem.csv
cat GeneMat.de.txt | grep -E 'Pattern2|Pattern4|Pattern5' > DEG_chrA-chrD_EBSeq_rsem.csv
cat GeneMat.de.txt | grep -E 'Pattern2|Pattern3|Pattern5' > DEG_chrB-chrD_EBSeq_rsem.csv


