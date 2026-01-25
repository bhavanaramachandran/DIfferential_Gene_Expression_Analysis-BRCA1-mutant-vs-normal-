# Differential Gene Expression Analysis: BRCA1-mutant ER negative vs Normal

## Objective
Identify differentially expressed genes and enriched pathways between **BRCA1-mutant ER negative tumors vs normal breast tissue** samples using bulk RNA-seq and **DESeq2**.

## Data
- Source: SRA Database (Shah et. al., 2022)
- Accession no. : PRJNA751555
- Samples: BRCA1-mutant (n = 9), Normal (n = 6)
- Notes on inclusion/exclusion: Analysis restricted to BRCA1 germline mutation primary tumor samples with ER negative status.

## Pipeline 
1. Input: counts matrix source: featureCounts
2. QC: Adapter trimming and removal of low quality reads with Trimmomatic
3. Differential expression: DESeq2  (Wald test; unshrunken results used for GSEA ranking; log2 fold-change shrinkage applied for visualization)
4. Annotation: GENCODE v19; GRCh37.p13.genome.fa
5. Visualization: Volcano plot and heatmap
6. Enrichment: GSEA with Hallmark geneset

## Key results 
- Main finding 1: The most significantly upregulated
Genes in BRCA1 mutant ER negative samples include NDC80 kinetochore complex component, Centromere Protein F, Methylenetetrahydrofolate dehydrogenase 2, and the most significantly downregulated genes were Peroxisome Proliferator-Activated Receptor Gamma, Cysteine_Rich Transmembrane BMP regulator 1, and Dermatopontin.

- Main finding 2:  Genes associated with
MTORC1 signaling, Myc, G2M Checkpoint, E2F targets and glycolysis are the most
significantly enriched as upregulated pathways and a downregulation of
genes associated with IL6/JAK/STAT3 signaling and inflammatory signaling was
observed.


- Main finding 3: Results broadly agree with known biological changes in BRCA1 ER negative tumors.

Figures:
- Volcano plot: `figures/volcano.png`
- PCA plot: `figures/pca.png`
- Enrichment: `figures/gsea_top.png`

## Reproducibility (minimal)
Note: This workflow reflects the exact commands and scripts used in this analysis and is intended as a transparent, learning-focused pipeline rather than a full production workflow.

To rerun:
1. Install R (≥ 4.1.2), Python(≥ 3.13.2)
2. Install packages listed in ‘Software/packages’ section below
3. The main scripts (excluding visualization) can be run in this order:
   - `scripts/prefetch_controls.sh and scripts/prefetch_samples.sh`
   - `scripts/fasterqdump_controls.sh and. scripts/fasterqdump_samples.sh`
   - `scripts/fastqc_pretrimming.sh`
   - `scripts/trimmomatic.sh`
   - `scripts/fastqc_posttrimming.sh`
   - `scripts/STAR_download.sh`
   - `scripts/index_ref.sh`
   - `scripts/alignment.sh`
   - `scripts/featurecounts.sh`
   - `scripts/analysis.R (DESeq2)`
   - `scripts/significant_results.R`
   - `scripts/run_gsea.py`

### Software / packages
- R version: 4.1.2
- Key packages: DESeq2 (v1.34.0),  STAR aligner (v2.7.10b), featureCounts (v2.0.3), fastQC(v0.11.9), Trimmomatic (v0.39), Gseapy (v1.1.8). Prefetch (v3.2.0), Fasterqdump (v3.2.0), Python(v3.13.2)

## Outputs
- Differential expression results: results/DESeq2_results_with_symbols.csv, results/significant_downregulated_genes.csv, results/significant_upregulated_genes.csv
- Figures: results/volcano_plot.png, results/heatmap_top_100_DE_genes.png, results/top_enriched_barplot.png
## Detailed analysis report
A comprehensive technical report describing the study design, quality control,
parameter choices, differential expression analysis, pathway enrichment, and
biological interpretation is available here:

- docs/Technical_Report_RNAseq_BRCA1_Reanalysis.pdf

**Note:** This report documents a learning-focused reanalysis of publicly available data (PRJNA751555) and is not a peer-reviewed manuscript or a novel biological claim.

##Acknowledgements:
I would like to thank Rajiv Gandhi Centre for Biotechnology, Trivandrum for providing
server access and computational support necessary for performing this analysis. 


