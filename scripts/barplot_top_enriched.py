import pandas as pd
import matplotlib.pyplot as plt

# === Load GSEA results ===
gsea_csv = "/Users/bhavanaramachandran/Desktop/gsea_results/hallmark/gseapy.gene_set.prerank.report.csv"
df = pd.read_csv(gsea_csv)

# === Drop NA and sort ===
df = df.dropna(subset=["NES"])
df_sorted = df.sort_values("NES", ascending=False)

# === Top 10 up and down pathways ===
top_up = df_sorted.head(10)
top_down = df_sorted.tail(10)

# === Combine and sort again for plotting ===
top_combined = pd.concat([top_up, top_down])
top_combined = top_combined.sort_values("NES")

# === Plot ===
plt.figure(figsize=(10, 6))
bars = plt.barh(top_combined["Term"], top_combined["NES"], color=["#d62728" if x > 0 else "#1f77b4" for x in top_combined["NES"]])
plt.xlabel("Normalized Enrichment Score (NES)")
plt.title("Top Enriched Gene Sets (Hallmark)")
plt.axvline(0, color='grey', linestyle='--')
plt.tight_layout()

# === Save plot ===
plt.savefig("/Users/bhavanaramachandran/Desktop/gsea_results/hallmark/top_enriched_barplot.png", dpi=300)
plt.show()

