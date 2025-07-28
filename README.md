# WES_Pipeline
# Whole Exome Sequencing (WES) Human Pipeline

## **Overview**
This repository contains a bioinformatics pipeline for analyzing Whole Exome Sequencing (WES) data of human samples. The pipeline is designed for accurate variant calling, annotation, and quality control using industry-standard tools such as **BWA-MEM, GATK4, Samtools, and FastQC**. 

The pipeline is built to handle **paired-end Illumina reads (2 × 100 bp or 150 bp)** and follows GATK Best Practices for germline variant discovery.

---

## **Pipeline Workflow**
The pipeline consists of the following steps:

1. **Quality Control (QC)**
   - Tool: FastQC / MultiQC
   - Purpose: Assess raw FASTQ read quality (e.g., adapter contamination, GC content).

2. **Read Trimming (Optional)**
   - Tool: Trimmomatic or fastp
   - Purpose: Remove low-quality bases and adapters.

3. **Alignment**
   - Tool: BWA-MEM
   - Input: Cleaned FASTQ reads
   - Output: Aligned BAM file.

4. **Post-Alignment Processing**
   - Sorting and indexing (Samtools)
   - Marking duplicates (Picard or GATK MarkDuplicates)
   - Base quality score recalibration (BQSR) using GATK.

5. **Variant Calling**
   - Tool: GATK HaplotypeCaller (GVCF mode)
   - Output: Raw variant calls (VCF files).

6. **Variant Filtering**
   - Apply GATK Variant Quality Score Recalibration (VQSR) or hard filtering.

7. **Variant Annotation**
   - Tool: Annovar, SnpEff, or VEP
   - Output: Annotated VCF with functional impact predictions.

8. **Coverage Analysis (Optional)**
   - Tool: mosdepth or bedtools
   - Generate per-target or per-gene coverage statistics.

---

## **Repository Contents**
- `wes_pipeline.sh` – Main pipeline script (bash).
- `requirements.txt` – List of software/tools needed.
- `config/` – Configuration files (e.g., reference genome, BED targets).
- `example_data/` – Example FASTQ or test dataset (if available).
- `README.md` – This documentation.

---

## **Requirements**
- **OS:** Linux (tested on Ubuntu/WSL)
- **Tools:**
  - BWA
  - Samtools
  - GATK4
  - Picard
  - FastQC, MultiQC
  - bcftools
  - Python 3 / R (for downstream analysis)
  
Reference genome: **hg38 (or hg19)** with .fa, .fai, and .dict files.

---

## **Usage**
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/WES_Pipeline.git
   cd WES_Pipeline
