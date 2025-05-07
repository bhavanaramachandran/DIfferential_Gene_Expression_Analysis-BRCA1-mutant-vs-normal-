# volcano_plot.R
# Create volcano plot of DESeq2 results

library(ggplot2)
library(dplyr)

# Load results
res <- read.csv("results/deseq2_results_with_symbols.csv")

# Basic filtering for plotting
res <- res %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE                              ~ "Not significant"
    )
  )

# Label top 20 most up and 20 most downregulated
top_genes <- res %>%
  filter(!is.na(gene_name)) %>%
  arrange(padj) %>%
  filter(significance != "Not significant") %>%
  slice_max(log2FoldChange, n = 20) %>%
  bind_rows(
    res %>%
      filter(!is.na(gene_name)) %>%
      arrange(padj) %>%
      filter(significance != "Not significant") %>%
      slice_min(log2FoldChange, n = 20)
  )

# Volcano plot
p <- ggplot(res, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  geom_text(data = top_genes, aes(label = gene_name), size = 3, vjust = 1.3, check_overlap = TRUE) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: BRCA1-mutant vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  )

# Save plot
ggsave("results/volcano_plot.png", plot = p, width = 10, height = 7, dpi = 300)
