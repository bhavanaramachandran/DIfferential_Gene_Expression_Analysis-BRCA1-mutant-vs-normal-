#!/bin/bash

# Move to the location of the STAR output files
cd "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/STAR_output" || exit

# Set paths
STARfiles="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/STAR_output"
outdir="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/FeatureCounts_output"

# Make output directory if it doesn't exist
mkdir -p "$outdir"

# Run featureCounts
featureCounts \
-a /home/bhavana/gencode.v19.annotation.gtf \
-o "$outdir/counts.txt" \
-T 8 \
-p --countReadPairs -B -C \
-t exon -g gene_id \
"$STARfiles"/*_Aligned.sortedByCoord.out.bam
