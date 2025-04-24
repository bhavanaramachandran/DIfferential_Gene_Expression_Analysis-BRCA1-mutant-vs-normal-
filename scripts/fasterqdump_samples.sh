#!/bin/bash

SAMPLES="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/prefetch_samples"

find "$SAMPLES" -type f -name "*.sra" | while read -r sra_file; do
    echo "Processing $sra_file"
    fasterq-dump "$sra_file"
done
