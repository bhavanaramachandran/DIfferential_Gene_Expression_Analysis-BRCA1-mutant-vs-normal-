# Differential Gene Expression Analysis: BRCA1-mutant vs Normal

This repository contains the RNA-seq pipeline to identify differentially expressed genes between BRCA1-mutant and normal samples.

## Tools
- SRA Toolkit
- Trimmomatic
- STAR
- featureCounts
- DESeq2

## Steps
1. Download raw reads from SRA
2. Quality control and trimming
3. Alignment to reference genome (GRCh37/hg19)
4. Quantification with featureCounts
5. Differential expression analysis with DESeq2
