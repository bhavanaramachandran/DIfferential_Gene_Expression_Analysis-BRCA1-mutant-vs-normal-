# significant_results.R

# Load libraries
library(dplyr)

# Read annotated DESeq2 results
res <- read.csv("results/deseq2_results_with_symbols.csv")

# Define significant thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Filter for significantly upregulated genes
up <- res %>%
  filter(padj < padj_cutoff, log2FoldChange > log2fc_cutoff)

# Filter for significantly downregulated genes
down <- res %>%
  filter(padj < padj_cutoff, log2FoldChange < -log2fc_cutoff)

# Save to CSV
write.csv(up, "results/significant_upregulated_genes.csv", row.names = FALSE)
write.csv(down, "results/significant_downregulated_genes.csv", row.names = FALSE)
