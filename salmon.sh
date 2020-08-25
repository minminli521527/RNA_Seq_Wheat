### 1 Salmon quantitative
# 1.1 Index the transcript
# The parameter-i indicates for index, followed by the index name. If you are still a species in the future, you can use the index established this time directly without repeating this step.
# Triticum
kallisto index /kallisto/0_data/Triticum_aestivum.IWGSC.cdna.all.fa.gz -i /kallisto/0_data/Triticum_aestivum.index

# 1.2 Identification and quantification of transcripts
# Paired-end sequencing
# -b: The number of samples for the self-service method, the author recommends to be set to at least 30, if the downstream does not need to use the booststraps value (for example, do not do sleuth, QTL analysis), -b can be set to 0 to speed up the operation.

# index
salmon index -p 4 -t all_merged.fasta -i all_merged_salmon_index
# Quantitative
for i in {CS_CT_1_1,CS_CT_2_1,CS_CT_3_1,ABA1h_1_1,ABA1h_2_1,ABA1h_3_1,ABA12h_1_1,ABA12h_2_1,ABA12h_3_1,ABA24h_1_1,ABA24h_2_1,ABA24h_3_1};do salmon quant -i Uall_merged_salmon_index -l A -1 ../clean/${i}_1_clean.fq.gz -2 ../clean/${i}_2_clean.fq.gz -p 8 --numBootstraps 100 -o ${i}_quant;done

# Parameter-- numBootstraps 100 must be put in. 
# put the generated * _ quant folder under the salmon_result folder

mkdir salmon_result
mv *_quant salmon_result/

# At this point, the quantitative part is completed.




conda install -c bioconda r-wasabi -y
R
library(wasabi)
sfdirs <- file.path("salmon_result", c("CS_CT_1_quant","CS_CT_2_quant","CS_CT_3_quant","ABA1h_1_quant","ABA1h_2_quant","ABA1h_3_quant","ABA12h_1_quant","ABA12h_2_quant","ABA12h_3_quant","ABA24h_1_quant","ABA24h_2_quant","ABA24h_3_quant"))
prepare_fish_for_sleuth(sfdirs)
