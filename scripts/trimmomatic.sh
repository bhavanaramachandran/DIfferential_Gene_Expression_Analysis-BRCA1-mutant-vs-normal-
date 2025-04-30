#!/bin/bash

# Set paths
fastq_dir="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/all_fastq_files/"
trimmed_dir="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/trimmed_fastq/"
trimmomatic_path="/home/bhavana/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_file="/home/bhavana/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# Create output directory if it doesn't exist
mkdir -p "$trimmed_dir"

# Loop through each _1.fastq file
for file in "$fastq_dir"/*_1.fastq; do
    base_name=$(basename "$file" _1.fastq)
    R1="$fastq_dir/${base_name}_1.fastq"
    R2="$fastq_dir/${base_name}_2.fastq"

    out_R1_paired="$trimmed_dir/${base_name}_1_paired_trimmed.fastq"
    out_R1_unpaired="$trimmed_dir/${base_name}_1_unpaired.fastq"
    out_R2_paired="$trimmed_dir/${base_name}_2_paired_trimmed.fastq"
    out_R2_unpaired="$trimmed_dir/${base_name}_2_unpaired.fastq"

    if [[ -f "$R2" ]]; then
        echo "Trimming $R1 and $R2..."
        java -jar "$trimmomatic_path" PE -phred33 "$R1" "$R2" \
            "$out_R1_paired" "$out_R1_unpaired" \
            "$out_R2_paired" "$out_R2_unpaired" \
            ILLUMINACLIP:"$adapter_file":2:20:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:12 MINLEN:30
    else
        echo "No matching pair for $R1, skipping."
    fi
done
