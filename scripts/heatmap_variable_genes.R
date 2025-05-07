# heatmap_variable_genes.R
# Heatmap of top 100 most variable genes (with gene symbols)

library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Load normalized counts and DESeq2 results
norm_counts <- read.csv("results/normalized_counts.csv", row.names = 1)
res_annotated <- read.csv("results/deseq2_results_with_symbols.csv")

# Strip version suffix from gene_id for matching
rownames(norm_counts) <- gsub("\\.\\d+$", "", rownames(norm_counts))
res_annotated$gene_id <- gsub("\\.\\d+$", "", res_annotated$gene_id)

# Add gene symbols as rownames where available
gene_symbols <- res_annotated %>%
  select(gene_id, gene_name) %>%
  filter(!is.na(gene_name)) %>%
  distinct()

# Join symbols to counts
norm_counts$gene_id <- rownames(norm_counts)
merged <- inner_join(norm_counts, gene_symbols, by = "gene_id")

# Remove versionless Ensembl ID column
rownames(merged) <- make.unique(merged$gene_name)
merged <- merged %>% select(-gene_id, -gene_name)

# Compute variance and select top 100 variable genes
merged$variance <- apply(merged, 1, var)
top_genes <- merged %>%
  arrange(desc(variance)) %>%
  slice(1:100) %>%
  select(-variance)

# Z-score scaling
scaled <- t(scale(t(as.matrix(top_genes))))

# Column annotation (mutant = 9, normal = 6)
sample_names <- colnames(scaled)
condition <- c(rep("mutant", 9), rep("normal", 6))
annotation_col <- data.frame(condition = condition)
rownames(annotation_col) <- sample_names

# Generate heatmap
pheatmap(scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         filename = "results/heatmap_top_100_variable_genes.png",
         width = 10,
         height = 12)
