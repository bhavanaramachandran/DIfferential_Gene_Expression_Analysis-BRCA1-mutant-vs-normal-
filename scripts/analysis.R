# analysis.R
# Author: Bhavana R
# Purpose: DE analysis of BRCA1-mutant vs normal with gene symbol annotation

# Load libraries
library(DESeq2)
library(dplyr)
library(stringr)
library(rtracklayer)  # for reliable GTF parsing

# ---- Step 1: Read featureCounts output ----
counts <- read.delim(
  "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/data/FeatureCounts_output/counts.txt",
  comment.char = "#", check.names = FALSE)

count_data <- counts[, 7:ncol(counts)]
colnames(count_data) <- gsub(".*/|_Aligned.sortedByCoord.out.bam", "", colnames(count_data))
rownames(count_data) <- counts$Geneid

# ---- Step 2: Define sample metadata ----
sample_names <- colnames(count_data)
condition <- c(rep("mutant", 9), rep("normal", 6))
col_data <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("normal", "mutant")))

# ---- Step 3: Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef = "condition_mutant_vs_normal", type = "apeglm")
res_ordered <- res[order(res$padj), ]

# ---- Step 4: Extract gene symbol mapping using rtracklayer ----
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

# ---- Step 5: Annotate DE results ----
res_ordered_df <- as.data.frame(res_ordered) %>%
  mutate(gene_id = rownames(.))

# Clean version suffix from DESeq2 gene_id for matching
res_ordered_df$gene_id <- gsub("\\.\\d+$", "", res_ordered_df$gene_id)

res_annotated <- left_join(res_ordered_df, gene_map, by = "gene_id") %>%
  select(gene_id, gene_name, everything())

# ---- Step 6: Save results ----
dir.create("/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/results", showWarnings = FALSE)

write.csv(
  res_annotated,
  file = "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/results/deseq2_results_with_symbols.csv",
  row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)

write.csv(
  as.data.frame(norm_counts),
  file = "/home/bhavana/DIfferential_Gene_Expression_Analysis-BRCA1-mutant-vs-normal-/results/normalized_counts.csv")
