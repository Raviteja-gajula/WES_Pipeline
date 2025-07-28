 #!/bin/bash
###################################################### VARIANT CALLING STEPS ####################################################################
# directories
ref="/home/Ravi2023/NGS_test/hg38.fa"
known_sites="/home/Ravi2023/NGS_test/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/Ravi2023/NGS_test/aligned_reads"
reads="/home/Ravi2023/NGS_test/reads"
trimmed_reads="/home/Ravi2023/NGS_test/reads/trimmed_reads"
results="/home/Ravi2023/NGS_test/results"
data="/home/Ravi2023/NGS_test/data"

# ------------------
# STEP 1: QC - Run fastqc 
# -------------------
echo "STEP 1: QC - Run fastqc"
fastqc ${reads}/Yoursample_R1_001.fastq.gz -o ${reads}/
fastqc ${reads}/Yoursample_R2_001.fastq.gz -o ${reads}/


# ------------------
# STEP 2: Trimming
# -------------------
echo "Step 2: Trimming the adapter"

R1="/home/Ravi2023/NGS_test/reads/PG250425111536_R1_001.fastq.gz"     
R2="/home/Ravi2023/NGS_test/reads/PG250425111536_R2_001.fastq.gz"      
output="/home/Ravi2023/NGS_test/reads/trimmed_reads"
CPU_CORES=4

adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


trimmed_R1="${output}/$(basename $R1 .fastq.gz)_trimmed.fastq.gz"
trimmed_R2="${output}/$(basename $R2 .fastq.gz)_trimmed.fastq.gz"
report_file="${output}/cutadapt_report.txt"

echo "Starting trimming process..."
echo "  Read 1: $R1"
echo "  Read 2: $R2"

cutadapt \
  -g $adapter1 \
  -G $adapter2 \
  -o $trimmed_R1 \
  -p $trimmed_R2 \
  -j $CPU_CORES \
  --minimum-length 20 \
  $R1 $R2 > $report_file

echo "Trimming complete!"
echo "Results saved in: $output"
echo "Trimmed files:"
echo "  $trimmed_R1"
echo "  $trimmed_R2"
echo "Report: $report_file"

# ------------------
# STEP 3: QC - Run fastqc 
# -------------------
echo "STEP 3: QC - Run fastqc"
fastqc ${trimmed_reads}/Yoursample_R1_001_trimmed.fastq.gz -o ${reads}/
fastqc ${trimmed_reads}/Yoursample_R2_001_trimmed.fastq.gz -o ${reads}/

# --------------------------------------
# STEP 4: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 4: Map to reference using BWA-MEM"

# BWA index reference 
#bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:Yoursample\tPL:ILLUMINA\tSM:Yoursample" ${ref} ${trimmed_reads}/Yoursample_R1_001_trimmed.fastq.gz ${trimmed_reads}/Yoursample_R2_001_trimmed.fastq.gz > ${aligned_reads}/PG250425111536.paired.sam


# -----------------------------------------
# STEP 5: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 5: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/Yoursample.paired.sam -O ${aligned_reads}/Yoursample_sorted_dedup_reads.bam

# ----------------------------------
# STEP 6: Base quality recalibration
# ----------------------------------


echo "STEP 6: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/Yoursample_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/Yoursample_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/Yoursample_sorted_dedup_bqsr_reads.bam 

# -----------------------------------------------
# STEP 7: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 7: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/Yoursample_sorted_dedup_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/Yoursample_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

# ----------------------------------------------
# STEP 8: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 8: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/Yoursample_sorted_dedup_bqsr_reads.bam -O ${results}/Yoursample_raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/Yoursample_raw_variants.vcf --select-type SNP -O ${results}/Yoursample_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/Yoursample_raw_variants.vcf --select-type INDEL -O ${results}/Yoursample_raw_indels.vcf


# ----------------------------------------
# STEP 9: Variant Filtering - SNPs & INDELs
# ----------------------------------------

echo "STEP 9: Variant Filtering - SNPs & INDELs"

# Filter SNPs
gatk VariantFiltration \
 -R ${ref} \
 -V ${results}/Yoursample_raw_snps.vcf \
 --filter-expression "QD < 2.0" --filter-name "QD2" \
 --filter-expression "FS > 60.0" --filter-name "FS60" \
 --filter-expression "MQ < 40.0" --filter-name "MQ40" \
 --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
 -O ${results}/Yoursample_filtered_snps.vcf

# Filter INDELs
gatk VariantFiltration \
 -R ${ref} \
 -V ${results}/Yoursample_raw_indels.vcf \
 --filter-expression "QD < 2.0" --filter-name "QD2" \
 --filter-expression "FS > 200.0" --filter-name "FS200" \
 --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
 -O ${results}/Yoursample_filtered_indels.vcf


# ----------------------------
# STEP 10: Merge VCFs
# ----------------------------

echo "STEP 10: Merge SNP & INDEL VCFs"

gatk MergeVcfs \
 -I ${results}/Yoursample_filtered_snps.vcf \
 -I ${results}/Yoursample_filtered_indels.vcf \
 -O ${results}/Yoursample_filtered_all.vcf


# ----------------------------------------
# STEP 11: Extract only PASS variants
# ----------------------------------------

echo "STEP 11: Extract only PASS variants"

# Ensure bcftools is installed and accessible
bcftools view -f PASS ${results}/Yoursample_filtered_all.vcf -o ${results}/Yoursample_final_PASS_variants.vcf

echo "Filtering complete. Final PASS variants saved at:"
echo "${results}/final_PASS_variants.vcf"

echo "WES pipeline completed successfully!"
