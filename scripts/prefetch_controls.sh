#!/bin/bash
SRA_LIST="/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/control_sra_ids.txt"  #stores the SRA list file into a variable
CONTROL="//home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/prefetch_controls"       #Ensure you already have made this folder

while read -r SRA_ID; do              #reads each line in the list into the variable SRA_ID
    SRA_DIR="$CONTROL/$SRA_ID"     # -r ensures that backslashes in SRA IDs are treated literal.
                                                                 #SRA_DIR creates unique subfolder for each SRA ID
    prefetch "$SRA_ID" -O "$SRA_DIR"               #downloads each SRA ID .sra file into output  
done < "$SRA_LIST" 
