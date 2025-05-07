# heatmap_de_genes.R
# Heatmap of top 100 most significantly differentially expressed genes

library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Load normalized counts and DE results
norm_counts <- read.csv("results/normalized_counts.csv", row.names = 1)
res <- read.csv("results/deseq2_results_with_symbols.csv")

# Clean gene IDs (remove version suffix)
rownames(norm_counts) <- gsub("\\.\\d+$", "", rownames(norm_counts))
res$gene_id <- gsub("\\.\\d+$", "", res$gene_id)

# Filter top 100 DE genes by adjusted p-value
top_de_genes <- res %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:100) %>%
  filter(!is.na(gene_name))

# Extract their normalized expression
counts_top <- norm_counts[rownames(norm_counts) %in% top_de_genes$gene_id, ]

# Attach gene names as rownames
counts_top$gene_id <- rownames(counts_top)
counts_named <- inner_join(counts_top, top_de_genes[, c("gene_id", "gene_name")], by = "gene_id")
rownames(counts_named) <- make.unique(counts_named$gene_name)
counts_named <- counts_named %>% select(-gene_id, -gene_name)

# Z-score transform
scaled <- t(scale(t(as.matrix(counts_named))))

# Annotate columns
sample_names <- colnames(scaled)
condition <- c(rep("mutant", 9), rep("normal", 6))
annotation_col <- data.frame(condition = condition)
rownames(annotation_col) <- sample_names

# Draw heatmap
pheatmap(scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         filename = "results/heatmap_top_100_DE_genes.png",
         width = 10,
         height = 12)
