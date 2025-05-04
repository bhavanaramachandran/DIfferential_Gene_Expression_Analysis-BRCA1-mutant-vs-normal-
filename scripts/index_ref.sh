#!/bin/bash

cd /home/bhavana || exit

# Create required directories
mkdir -p STAR_index_files

# Run STAR index generation
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /home/bhavana/STAR_index_files \
     --genomeFastaFiles /home/bhavana/GRCh37.p13.genome.fa \
     --sjdbGTFfile /home/bhavana/gencode.v19.annotation.gtf \
     --sjdbOverhang 149 \
