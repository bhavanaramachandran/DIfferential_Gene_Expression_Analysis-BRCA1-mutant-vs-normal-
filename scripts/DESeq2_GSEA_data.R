# analysis_unshrunken.R
# Author: Bhavana R
# Purpose: DE analysis of BRCA1-mutant vs normal with gene symbol annotation for GSEA (UNSHRUNKEN RESULTS ONLY)

# ---- Load libraries ----
library(DESeq2)
library(dplyr)
library(stringr)
library(rtracklayer)  # for reliable GTF parsing
library(tibble)

# ---- Step 1: Read featureCounts output ----
counts <- read.delim(
  "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/FeatureCounts_output/counts.txt",
  comment.char = "#",
  check.names = FALSE
)

count_data <- counts[, 7:ncol(counts)]
colnames(count_data) <- gsub(".*/|_Aligned.sortedByCoord.out.bam", "", colnames(count_data))
rownames(count_data) <- counts$Geneid

# ---- Step 2: Define sample metadata ----
sample_names <- colnames(count_data)
condition <- c(rep("mutant", 9), rep("normal", 6))
col_data <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("normal", "mutant")))

# ---- Step 3: Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)

# ---- Step 4: Extract unshrunken DE results (includes 'stat' column) ----
res <- results(dds)

# ---- Step 5: Extract gene symbol mapping from GENCODE v19 GTF ----
extract_gene_map_rtracklayer <- function(gtf_file) {
  gtf <- rtracklayer::import(gtf_file)
  gene_df <- as.data.frame(gtf[gtf$type == "gene", ])

  tibble(
    gene_id = gsub("\\.\\d+$", "", gene_df$gene_id),
    gene_name = gene_df$gene_name
  ) %>%
    distinct() %>%
    filter(!is.na(gene_id) & !is.na(gene_name))
}

gene_map <- extract_gene_map_rtracklayer("/home/bhavana/gencode.v19.annotation.gtf")

# ---- Step 6: Annotate DE results ----
res_df <- as.data.frame(res) %>%
  mutate(gene_id = gsub("\\.\\d+$", "", rownames(.))) %>%
  left_join(gene_map, by = "gene_id") %>%
  select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# ---- Step 7: Save results for GSEA ----
write.csv(
  res_df,
  file = "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/GSEA_DESeq2_results.csv",
  row.names = FALSE
)

message("✅ Unshrunken DESeq2 results with stat column saved for GSEA.")
