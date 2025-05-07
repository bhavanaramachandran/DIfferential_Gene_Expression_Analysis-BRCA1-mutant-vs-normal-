#!/bin/bash

# Set paths
fastq_dir="/home/bhavana/trimmed_fastq"
output_dir="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/STAR_output"
genome_dir="/home/bhavana/STAR_index_files"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each paired FASTQ file (1st pair)
for file1 in "$fastq_dir"/*_1_paired_trimmed.fastq; do

    # Derive base sample name (removes the _1_paired_trimmed.fastq)
    base_name=$(basename "$file1" _1_paired_trimmed.fastq)

    # Construct path to matching second FASTQ file
    file2="$fastq_dir/${base_name}_2_paired_trimmed.fastq"

    # Set unique output prefix for each sample
    out_prefix="$output_dir/${base_name}_"

    # Run STAR alignment
    STAR --genomeDir "$genome_dir" \
         --readFilesIn "$file1" "$file2" \
         --runThreadN 8 \
         --outFileNamePrefix "$out_prefix" \
         --outSAMtype BAM SortedByCoordinate
done
