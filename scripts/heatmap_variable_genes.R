# heatmap_variable_genes.R
# Heatmap of top 100 most variable genes based on normalized expression

library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Load normalized counts
norm_counts <- read.csv("results/normalized_counts.csv", row.names = 1)

# Step 1: Compute variance across samples for each gene
norm_counts$variance <- apply(norm_counts, 1, var)

# Step 2: Select top 100 most variable genes
top_genes <- norm_counts %>%
  arrange(desc(variance)) %>%
  slice(1:100) %>%
  select(-variance)  # remove helper column

# Step 3: Z-score transformation by gene (row-wise scaling)
scaled <- t(scale(t(as.matrix(top_genes))))

# Step 4: Create sample annotation
sample_names <- colnames(scaled)
condition <- c(rep("mutant", 9), rep("normal", 6))
annotation_col <- data.frame(condition = condition)
rownames(annotation_col) <- sample_names

# Step 5: Generate heatmap
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
