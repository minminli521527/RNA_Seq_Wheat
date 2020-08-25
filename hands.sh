
### 1 Software installation.
### sratoolkit
curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz
tar xvf sratoolkit.2.8.2-ubuntu64.tar.gz
cd sratoolkit.2.8.2-ubuntu64/bin      
echo "export PATH=$PATH:/software/sratoolkit.2.8.2-ubuntu64/bin" >> ~/.bashrc   #Add environment variables.
source ~/.bashrc
### BWA
# https://sourceforge.net/projects/bio-bwa/files/
tar jxvf bwa-0.5.9rc1.tar.bz2
cd bwa-0.7.17/
make    # compilation
echo "export PATH=/software/bwa-0.7.17:$PATH" >> ~/.bashrc   #Add environment variables.
source ~/.bashrc
### GSNAP-GMAP
# http://research-pub.gene.com/gmap/src/README
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz
tar xf gmap-gsnap-2018-07-04.tar.gz
cd gmap-2018-07-04/
./configure --prefix=$HOME/opt/biosoft/gmap
make -j 20
make install
### BamBam
# https://sourceforge.net/p/bambam/activity/
# Read carefully the install and README in the installation package.
# First install SAMtools (including HTSlib), BAMtools, and Zlib.
# Add a LIBRARY_PATH.

echo "export LIBRARY_PATH=/programs/samtools/samtools-1.9:/software/RSEM-1.3.3/samtools-1.3/htslib-1.3:/software/eagle/htslib:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/rnaseq/lib:/software/zlib-1.2.11:/programs/samtools/samtools-1.9/win32:/programs/velvet/v1.2.10/x86_64/third-party/zlib-1.2.3:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/rnaseq/lib:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/software/bamtools/build/src/api:$LIBRARY_PATH" >> ~/.bashrc
echo "export LIBRARY_PATH=/software/bamtools/build/src/utils:$LIBRARY_PATH" >> ~/.bashrc
source ~/.bashrc
### seqtk
git clone https://github.com/lh3/seqtk.git 
cd seqtk
make
echo "export PATH=$PATH:/software/seqtk" >> ~/.bashrc   #Add environment variables.
source ~/.bashrc



### 2 Prepare the data to be analyzed
### Download
#  HANDS: a tool for genome-wide discovery of subgenome- specific base-identity in polyploids. BMC Genomics 14: 653. doi:10.1186/1471-2164-14-653.
# [Download](https://www.ncbi.nlm.nih.gov/sra/?term=SRA097144&utm_source=gquery&utm_medium=search)
# SRR949820.sra
# SRR949821.sra
# SRR949822.sra
# SRR949823.sra
### Extract
# --split-3: If it is two-end sequencing data, two files are output, and if not only one file is output.
# --gzip: Compressed files in gzip output format (fastqc software can directly identify gzip compressed files)
# -O-outdir: Output to the specified folder.
# -A: File name after decompression.
fastq-dump/HANDS/SRR949820/SRR949820.sra--gzip--split-3-O. /raw_data.
### Extract some of the data for analysis.



head -n 400000 SRR949820_1.fastq.gz > SRR949820.1_1.fastq.gz
head -n 400000 SRR949820_2.fastq.gz > SRR949820.1_2.fastq.gz
head -n 400000 SRR949821_1.fastq.gz > SRR949821.1_1.fastq.gz
head -n 400000 SRR949821_2.fastq.gz > SRR949821.1_2.fastq.gz
head -n 400000 SRR949822_1.fastq.gz > SRR949822.1_1.fastq.gz
head -n 400000 SRR949822_2.fastq.gz > SRR949822.1_2.fastq.gz
head -n 400000 SRR949823_1.fastq.gz > SRR949823.1_1.fastq.gz
head -n 400000 SRR949823_2.fastq.gz > SRR949823.1_2.fastq.gz



### 3 Data quality filtering
### (1)fastqc
# Options:
# -o --outdir To output the directory, you need to create your own directory.
# -noextract: Do not extract files.
# -f: Specifies the type of input file, supports fastq|bam|sam files in three formats, and is automatically recognized by default.

ls *.fastq.gz | while read id;
do
fastqc -o /1_clean_reads/fastqc_1/ -t 3 -noextract -f fastq $id
done
### (2) Excision joint sequence
# Options:
# PE: Filter paired-end sequencing data (if the data is single-end sequencing, use SE).
# phred33: The quality value format of the Fastq file is phred33. Generally, the format of the second-generation sequencing data is basically phred33. If you are not sure, you can consult the sequencing company.
# threads: Set the number of threads.
# seq*.fq.gz: Fastq files to be filtered.
# seq*.clean.fq.gz: Fastq file after filtering.
# ILLUMINACLIP: ./adapters/TruSeq3-PE.fa refers to the removal of the TruSeq3 adapter sequence under the Illumina sequencing platform. For specific adapter sequences, please consult the sequencing company. The 3 numbers (2:30:10) following the adaptor sequence indicate that two base mismatches are allowed when aligning the adaptor sequence. If the match rate between the two reads of paired-end sequencing and the adaptor sequence exceeds 30% , It will be cut off, and if the match rate of a single read with the linker sequence exceeds 10%, it will also be cut off.
# SLIDINGWINDOW:4:15: indicates that the 4bp window is used for sliding window statistics, and the sequence of the window whose average quality of excision bases is lower than 15 and beyond.
# LEADING:3: Remove leading low quality or N bases (lower quality 3).
# TRAILING:3: Remove the tailing low quality or N bases (lower than quality 3).
# MINLEN:36: means to remove reads with a length less than 36 after filtering.

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
fastqc -o /software/HANDS/20200525/1_clean_reads/fastqc_2 -t 3 -noextract -f fastq $id
done



### 4 Use B subgenomic cdna to construct reference sequence
### (4.1) Extract wheat B genome cdna data
# Prepare to download wheat cdna data Triticum_aestivum.IWGSC.cdna.all.fa
# Separate by " ", display the first column (remove the characters after each gene name)
cut -d " " -f 1 Triticum_aestivum.IWGSC.cdna.all.fa > Triticum_aestivum.IWGSC.cdna1.all.fa

# Extract the transcript ID of the cdna sequence B genome
less Triticum_aestivum.IWGSC.cdna1.all.fa | grep TraesCS.B | sed 's/>//g' > Triticum_aestivum.IWGSC.cdna.B_ID

# Extract sequence based on gene ID list
# Method 1: Extract with seqtk subseq tool under linux system
seqtk subseq Triticum_aestivum.IWGSC.cdna1.all.fa Triticum_aestivum.IWGSC.cdna.B_ID > Triticum_aestivum.IWGSC.cdna.B.fa1
# Method 2: Use python script: 1.py: The content is as follows:

python3 ../Script/1.py
*************************************************************************************************************************************
### (1) Extract sequence ID name (de-duplication)
gene_list = [] # Set an empty list for storing protein names
with open(r'/software/HANDS/20200627/0_data/Triticum_aestivum.IWGSC.cdna.B_ID','r') as a:  #Read the file containing the protein name
    for i in a: # for loop to traverse line by line
        if not i.startswith('#'): # Skip comment lines in the file
            gene_id = i.split(' ')[0] # Press the space to separate each line of information and select the first column of protein name
            if gene_id not in gene_list: # If there is no current element in the temporary list, then append
                gene_list.append(gene_id) # append function is used to add qualified protein names to the list
print(len(gene_list)) # View the number of proteins

### (2) Extract the sequence of the corresponding protein from the fa file
final_seq = {} # Set an empty dictionary to store the results
parse_check = False # Set the variable to determine whether the read line is the line of the desired protein name, first set it to flase
with open(r'/software/HANDS/20200627/0_data/Triticum_aestivum.IWGSC.cdna1.all.fa','r') as a:
    for i in a:
        if i.startswith('>'): # Determine whether it starts with >
            seq_ID = i.split(' ')[0][1:] # Read the protein name of the fa file
            if seq_ID in gene_list: # Determine whether the protein name in the fa file is in the required protein list
                parse_check = True # If included in gene_list, parse_check will be re-assigned to True
                seq=[] # Set an empty list to store the sequence
                final_seq[seq_ID] = seq # Set the value corresponding to the key in the dictionary as a list
            else:
                parse_check = False # If the protein of fa is not in the gene_list, parse is still false
        elif parse_check == True: # If parse is true (ie fa protein is in gene_list)
            seq_part= i.split('\n')[0] # Add sequence to sequence list seq
            seq.append(seq_part)

# #Output the results in the dictionary to results.fa
outfile = open('Triticum_aestivum.IWGSC.cdna.B.fa', 'w')

for id in final_seq:
    outfile.write('>%s\n' % id) # Set string format to extract id information (ie protein name)
    sequence = final_seq[id] # Get the sequence corresponding to id
    for i in sequence:
        outfile.write('%s' % (i)) # Get sequence information and format output
    outfile.write('\n')
outfile.close()
*************************************************************************************************************************************
# Delete the blank line in the Triticum_aestivum.IWGSC.cdna.B.fa file
sed '/^\s*$/d' Triticum_aestivum.IWGSC.cdna.B.fa > Triticum_aestivum.IWGSC.cdna.B.fa2



### (4.2) Extract the longest transcript from wheat B genome cdna
### (4.2.1) Extract the ID name of the longest transcript sequence (get the longest_id.txt file)
# Use python script: 2.py: The content is as follows:

python3 ../Script/2.py
*************************************************************************************************************************************
#!/usr/bin/env python
# -*- coding=utf-8 -*-
 
'''
Extract the longest transcript ID in the sequence file
Need to modify the ###### position parameter and open directory
'''
 
import sys
import re
 
Fasta=open("/software/HANDS/20200627/0_data/Triticum_aestivum.IWGSC.cdna.B.fa","r")
Sequence={}
 
## (1) Create a dictionary of IDs and sequences
for line in Fasta.readlines():
  content=line.strip()
  if content.startswith(">"):
    nameS=content[1:]
    name=re.sub(" .*$","",nameS)       ###### Modify the header of the sequence.
    Sequence[name]=''
  else:
    Sequence[name]+=content
Fasta.close()
#print(Sequence.keys())
 
## (2) Extract the sequence length of each transcript to form three columns.
Out=open("/software/HANDS/20200627/0_data/gene_id_len.txt","w")
for i in Sequence.keys():
    #print(i+"\t"+str(len(Sequence[i]))+"\t"+i.split(".")[0]+"\n")
    Out.write(i+"\t"+str(len(Sequence[i]))+"\t"+i.split(".")[0]+"\n")   # The first column is transcript ID, the second column is transcript length, and the third column is gene ID.
Out.close()
 
## (3) This step prepares to sort the out files and find the longest transcript of each gene
from operator import itemgetter
gene_id_len=open("/software/HANDS/20200627/0_data/gene_id_len.txt","r")
table = []
for line in gene_id_len:
    col = line.strip().split("\t")
    col[0] = str(col[0])
    col[1] = int(col[1])
    col[2] = str(col[2])
    table.append(col)                                              # Build nested lists
#print(table)
table_sorted = sorted(table, key=itemgetter(2, 1),reverse=True)                 # Sorted by column index 2, 1.
#print(table_sorted)
output_file = open("/software/HANDS/20200627/0_data/gene_id_len_sorted.txt","w")
for row in table_sorted:                                           # Traverse a nested list after reading sorting
    row = [str(x) for x in row]                                    # The list needs to be converted to a string before it can be written to text.
    output_file.write("\t".join(row) + '\n')
output_file.close()
 
### (4) Look for the longest transcript.
input_file=open("/software/HANDS/20200627/0_data/gene_id_len_sorted.txt","r")
dict2={}
for line in input_file.readlines():
    col = line.strip().split("\t")
    col[0] = str(col[0])
    col[1] = int(col[1])
    col[2] = str(col[2])
    if col[2] not in dict2:
        dict2[col[2]]=col[0]
    else:
        dict2[col[2]]+= '\t'+col[0]
#print(dict2)                     
#print(dict2.values())
list_values=list(dict2.values())
#print(list_values)
result_file = open("/software/HANDS/20200627/0_data/longest_id.txt","w")
for line in list_values:
    col = line.strip().split("\t")
    col[0] = str(col[0])
    #print(col[0])
    result_file.write(col[0]+"\n")
result_file.close()
*************************************************************************************************************************************

### (4.2.2) Extract the longest transcript sequence
# Use seqtk subseq tool to extract the longest transcript sequence under linux system

seqtk subseq Triticum_aestivum.IWGSC.cdna.B.fa1 longest_id.txt > Triticum_aestivum.IWGSC.cdna.B.longest.fa



### (4.3) The hands2.jar seq2ref command will generate Ta.reference.fa and Ta.reference.fa.gff files.
# A diagrammatic representation of the T. aestivum reference sequence used in our analyses. The reference was created by concatenating 56,954 T. aestivum UniGene sequences (T. aestivum UniGene Build 60) such that two consecutive UniGenes were separated by a gap of 200 N’s. The reference was a total of 61,117,745 bases in length out of which 11,391,000 bases were ‘N’s used as separators.
# Options:
# -i <str>     :   Input sequence file (multifasta format)
# -o <str>     :   Output file
# -n <str>     :   Header for the in silico reference
# -g <int>     :   Gap size between two sequences (Default: 200)
java -jar /hands2/hands2.jar seq2ref -i Triticum_aestivum.IWGSC.cdna.B.longest.fa -o Ta.reference.fa -n Ta.reference



### 5 Generate vcf files by alignment with subgenomic sequences. 
### (1) BWA build index. 
# After the index is established, five files will be generated with the suffix bwt pac ann amb sa
bwa index Ta.reference.fa

### (2) BWA alignment
ls *_1.clean.fastq.gz |while read id
do
bwa mem ./Ta.reference.fa \
./$id ./${id%_*}_2.clean.fastq.gz > ${id%_*}.sam
done
### (3) SAM filtering
# The official website script: https://genomics.lums.edu.pk/software/hands2/hands2v1.1.1.tar.gz
# <path>/filter sam -i <sam_file> [other options]
# options
# -o filtered SAM file for output
# -g gff file
# -x max insert size (use -1 to select reads which map on same contig only; +ve value for reads having insert size less than the max_insert_size)
# -q mapping quality threshold
# -p paired data? TRUE/FALSE
# -b require both reads in a pair to be mapped TRUE/FALSE
ls *1.sam |while read id
do
/software/HANDS/20200525/filter_sam_vcf/filter sam -i $id -o ${id%.*}.filter.sam -g ./Ta.reference.fa.gff -q 20 -p TRUE -b TRUE
done
### (3) sam--bam
# sam--bam
ls *.filter.sam |while read id
do
samtools view -bS $id > ${id%.*}.bam
done

### (4) bam--sort
samtools sort SRR949820.1.filter.bam SRR949820.1.filter.sorted
samtools sort SRR949821.1.filter.bam SRR949821.1.filter.sorted
samtools sort SRR949822.1.filter.bam SRR949822.1.filter.sorted
samtools sort SRR949823.1.filter.bam SRR949823.1.filter.sorted

### (5) bam--index
samtools index SRR949820.1.filter.sorted.bam
samtools index SRR949821.1.filter.sorted.bam
samtools index SRR949822.1.filter.sorted.bam
samtools index SRR949823.1.filter.sorted.bam


### (6) vcf
# -f: to input the fasta reference sequence with index file;
# -g: output to bcf format;

samtools mpileup -guSDf ./Ta.reference.fa SRR949820.1.filter.sorted.bam | bcftools view -cvNg - > SRR949820.vcf
samtools mpileup -guSDf ./Ta.reference.fa SRR949821.1.filter.sorted.bam | bcftools view -cvNg - > SRR949821.vcf
samtools mpileup -guSDf ./Ta.reference.fa SRR949822.1.filter.sorted.bam | bcftools view -cvNg - > SRR949822.vcf
samtools mpileup -guSDf ./Ta.reference.fa SRR949823.1.filter.sorted.bam | bcftools view -cvNg - > SRR949823.vcf


### (7) vcf filter
### (7.1) samtools vcfutils.pl--varFilter conmmand 
ls *.vcf |while read id
do
/software/samtools/v0.1.20/x86_64/bin/vcfutils.pl varFilter -1 0 -4 0 -d 3 -D 50000 $id > ${id%.*}.filter_samtools.vcf
done
### (7.2) hands
# Official website script: https://genomics.lums.edu.pk/software/hands2/hands2v1.1.1.tar.gz
# <path>/filter vcf -i <vcf_file> [other options]
# -o filtered VCF file for output
# -q Variant quality threshold
ls *.filter_samtools.vcf |while read id
do
/software/HANDS/20200525/filter_sam_vcf/filter vcf -i $id -o ${id%.*}_hands.vcf -q 20
done


### (8) Calculate base coverage
# java -jar hands2.jar coverage -i example/brassica/napus.sam -o example/brassica/coverage.bc 
ls *.filter.sam |while read id
do
java -jar /software/hands2/v1.1.1/hands2/hands2.jar coverage -i $id -o ${id%.*}.coverage.bc
done


### (9) Assign homoeallelic base-identities
# java -jar hands2.jar assign -i example/wheat/polyploid.sam -g example/wheat/reference.fa.gff -hsp example/wheat/polyploid.hsp -snp1 example/wheat/diploid1.snp -snp2 example/wheat/diploid2.snp -snp3 example/wheat/diploid3.snp -bc example/wheat/polyploid.bc -bc1 example/wheat/diploid1.bc -bc2 example/wheat/diploid2.bc -bc3 example/wheat/diploid3.bc -out1 example/wheat/out1.vcf -out2 example/wheat/out2.vcf -out3 example/wheat/out3.vcf -m true 

java -jar /software/hands2/v1.1.1/hands2/hands2.jar assign -i ./SRR949820.1.filter.sam -g ./Ta.reference.fa.gff -hsp ./SRR949820.filter_samtools_hands.vcf -snp1 ./SRR949821.filter_samtools_hands.vcf -snp2 ./SRR949822.filter_samtools_hands.vcf -snp3 ./SRR949823.filter_samtools_hands.vcf -bc ./SRR949820.1.filter.coverage.bc -bc1 ./SRR949821.1.filter.coverage.bc -bc2 ./SRR949822.1.filter.coverage.bc -bc3 ./SRR949823.1.filter.coverage.bc -out1 SRR949821.subgenome.vcf -out2 SRR949822.subgenome.vcf -out3 SRR949823.subgenome.vcf -m true 



### 6 Build SNP_index index
### (6.1) GSNAP builds an index for genome files
mkdir genome
# -d represents the name of the database
# gmap_build -d genome reference.fa -D ./genome
# Genome files do not support compressed files

/software/gsnap/v20200408/gmap-2020-04-08/util/gmap_build -d genome Ta.reference.fa -D ./genome

### (6.2) Build snp-index
# (1) Combine vcf files in order by column
paste SRR949821.subgenome.vcf SRR949822.subgenome.vcf SRR949823.subgenome.vcf > subgenome.all.vcf
# (2) Use the python script 3.py to build the snp.txt file as follows:
paste SRR949821.subgenome.vcf SRR949822.subgenome.vcf SRR949823.subgenome.vcf > subgenome.all.vcf

python3 ../Script/3.py
*************************************************************************************************************************************
#!/usr/bin/env python3

str1 = []
fo = open("/software/HANDS/20200627/0_data/snp.txt", "w")
with open("/software/HANDS/20200627/0_data/subgenome.all.vcf", 'r') as file:
    for line in file.readlines():
#      print(line)
       if line[0] == "#":
           continue
       else:
#          print(line)
           str1 = line.strip('\n')
#          print (str1)
           list1 = str1.split('\t')
#          print (list1)
#          print (list1[4])
           if "." in list1[4]:
               list1[4] = list1[3]
           else:
               list1[4] = list1[4]
               
           if "." in list1 [12]:
               list1[12] = list1[11]
           else:
               list1[12] = list1[12]
               
           if "." in list1[20]:
               list1[20] = list1[19]
           else:
               list1[20] = list1[20]

           i=(list1[0],list1[1],list1[4],list1[12],list1[20])
           print ("\t".join(i) + "\n")
           fo.write("\t".join(i) + "\n")
fo.close()
*************************************************************************************************************************************
# (3) Add column names in the first line of snp.txt
 #Chr\tPos\tSRR949821.sorted.bam\tSRR949822.sorted.bam\tSRR949823.sorted.bam
sed '1i #Chr\tPos\tSRR949821.sorted.bam\tSRR949822.sorted.bam\tSRR949823.sorted.bam' snp.txt > snp1.txt

# (4) Build snp-index
/software/bambam-1.4/bambam/scripts/snp2gsnap.pl snp1.txt Ta.reference.fa | /software/gsnap/v20200408/gmap-2020-04-08/src/iit_store -o new_snps.iit
mkdir snp-index
/software/gsnap/v20200408/gmap-2020-04-08/src/snpindex -D ./genome -d genome -V ./snp-index -v new_snps new_snps.iit



### 7 GSNAP SNP-tolerant alignment
# Align the hexaploid reads to the B subgenome
/software/gsnap/v20200408/gmap-2020-04-08/src/gsnap --gunzip -n 1 -N 1 -Q -t 10 --merge-distant-samechr --use-sarray=0 -D ./genome -d genome -V ./snp-index -v new_snps -A sam SRR949820.1_1.clean.fastq.gz SRR949820.1_2.clean.fastq.gz > SRR949820.1.SNP-tolerant.sam

# sam-bam-sort-index (1.3.1 version samtools-conda)
samtools view -Sb -@ 10 SRR949820.1.SNP-tolerant.sam -o SRR949820.1.SNP-tolerant.bam
samtools sort SRR949820.1.SNP-tolerant.bam > SRR949820.1.SNP-tolerant.sorted.bam
samtools index SRR949820.1.SNP-tolerant.sorted.bam



### 6 PolyCat alignment
### (6.1) Using the previous snp-index file, polyCat is then divided into subgenomic bam.
# PolyCat generates three sub-genome bam files from the hex body bam file.
# -p 1 for pairend
/software/bambam-1.4/bambam/bin/polyCat -x 1 -p 1 -s ./snp-index SRR949820.1.SNP-tolerant.sorted.bam
# or
/software/bambam-1.4/bambam/bin/polyCat -x 1 -p 1 -s ./snp1.txt SRR949820.1.SNP-tolerant.sorted.bam

### (6.2) The snp-index file is generated using interSnp and then polyCat is divided into sub-genome bam.
# The snp-index file is generated by interSnp.
/software/bambam-1.4/bambam/bin/interSnp SRR949821.1.filter.sorted.bam SRR949822.1.filter.sorted.bam SRR949823.1.filter.sorted.bam 1>interSnp_snpindex
# PolyCat generates three sub-genome bam files from the hex body bam file.
# -p 1 for pairend
/software/bambam-1.4/bambam/bin/polyCat -x 1 -p 1 -s ./interSnp_snpindex SRR949820.1.SNP-tolerant.sorted.bam

### (6.3) igv visualization
### Import the genome
# Build index:samtools faidx genome.fa, get the index file called genome.fai.
# Upload the genome, the software does not support compression format: Genomes->Load Genome from File -> select genome.fa and open it.
### Open the BAM file obtained by mapping.
# sort BAM file：samtools sort sample.bam -o sample.sort.bam
# Build index fo BAM file：samtools index sample.sort.bam
# The BAM file after uploading the sort: File-Load from File-select sample.sort.bam and open it.
# See the results of the presentation.



### 7 Subgenomic reads count counts
# GFF file conversion GTF file (with python scriptgff_to_gtf.py, as follows:)

python2 gff_to_gtf.py Ta.reference.fa.gff > Ta.reference.fa.gtf
*************************************************************************************************************************************
#!/usr/bin/env python

#http://blog.nextgenetics.net/?e=27

#python myScript.py myFile.gff > myFile.gtf

import sys

inFile = open(sys.argv[1],'r')

for line in inFile:
  #skip comment lines that start with the '#' character
  if line[0] != '#':
    #split line into columns by tab
    data = line.strip().split('\t')

    ID = ''

    #if the feature is a gene 
    if data[2] == "gene":
      #get the id
      ID = data[-1].split('ID=')[-1].split(';')[0]

    #if the feature is anything else
    else:
      # get the parent as the ID
      ID = data[-1].split('Parent=')[-1].split(';')[0]
    
    #modify the last column
    data[-1] = 'gene_id "' + ID + '"; transcript_id "' + ID

    #print out this new GTF line
    print '\t'.join(data)
*************************************************************************************************************************************

# Subgenomic reads calculate counts.
htseq-count -f bam --stranded=no -r name SRR949820.1.SNP-tolerant.sorted.SRR949821.1.filter.sorted.bam Ta.reference.fa.gtf > SRR949820.1.SNP-tolerant.sorted.A.txt
htseq-count -f bam --stranded=no -r name SRR949820.1.SNP-tolerant.sorted.SRR949822.1.filter.sorted.bam Ta.reference.fa.gtf > SRR949820.1.SNP-tolerant.sorted.B.txt
htseq-count -f bam --stranded=no -r name SRR949820.1.SNP-tolerant.sorted.SRR949823.1.filter.sorted.bam Ta.reference.fa.gtf > SRR949820.1.SNP-tolerant.sorted.D.txt




