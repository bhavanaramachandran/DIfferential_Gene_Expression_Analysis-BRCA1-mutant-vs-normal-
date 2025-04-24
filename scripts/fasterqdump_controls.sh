#!/bin/bash

CONTROLS="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/prefetch_controls"

find "$CONTROLS" -type f -name "*.sra" | while read -r sra_file; do
    echo "Processing $sra_file"
    fasterq-dump "$sra_file"
done
