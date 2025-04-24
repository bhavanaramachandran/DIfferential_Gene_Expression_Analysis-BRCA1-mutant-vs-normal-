#!/bin/bash

# Run FastQC on all FASTQ files before trimming
# Output will go to fastqc_pretrimming directory

INPUT_DIR="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/all_fastq_files"
OUTPUT_DIR="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/fastqc_pretrimming"

mkdir -p $OUTPUT_DIR

fastqc -o $OUTPUT_DIR $INPUT_DIR/*.fastq
