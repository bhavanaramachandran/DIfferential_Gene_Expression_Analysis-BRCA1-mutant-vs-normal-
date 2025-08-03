# Differential Gene Expression Analysis: BRCA1-Mutant vs Normal

This repository contains an RNA-seq analysis pipeline to identify differentially expressed genes between BRCA1-mutant ER neagtive and normal samples using publicly available data.

## 📂 Data Source

The RNA-seq data used in this analysis was obtained from the NCBI Sequence Read Archive (SRA) under BioProject accession: **[PRJNA751555](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA751555)**.

This dataset includes RNA-seq samples from BRCA1-mutant ER negative and normal breast tissue. Raw FASTQ files were downloaded using the SRA Toolkit.




---

## Tools Used

- **SRA Toolkit** – for downloading raw sequencing data  
- **Trimmomatic** – for read trimming and quality control  
- **STAR** – for aligning reads to the human reference genome (GRCh37/hg19)  
- **featureCounts** – for gene-level quantification  
- **DESeq2** – for differential expression analysis in R
- **GSEApy** - for gene set enrichment analysis in python

---

##  Workflow Overview

1. **Download raw RNA-seq reads** using SRA Toolkit  
2. **Trim and clean reads** with Trimmomatic  
3. **Align reads** to GRCh37/hg19 using STAR  
4. **Quantify reads** with featureCounts  
5. **Analyze differential expression** with DESeq2
6. **Geneset enrichment analysis** with GSEApy

---

## Repository Structure

```text
├── data/          # Raw and trimmed FASTQ files, reference genome files  
├── scripts/       # Shell, R and Python scripts for all processing and analysis steps  
├── results/       # Output files: count matrix, DESeq2 results, plots  
└── README.md

