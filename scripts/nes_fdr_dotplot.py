import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Load GSEA result file ===
result_file = "/Users/bhavanaramachandran/Desktop/gsea_results/hallmark/gseapy.gene_set.prerank.report.csv"
df = pd.read_csv(result_file, sep=",")  # GSEA reports are usually tab-delimited

# === Filter: keep only terms with FDR < 0.25 ===
df = df[df['FDR q-val'] < 0.25]

# === Sort by NES ===
df = df.sort_values("NES", ascending=False)

# === Plot NES vs -log10(FDR) ===
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    x=df["NES"],
    y=-np.log10(df["FDR q-val"] + 1e-8),  # avoid log(0)
    c=df["NES"],
    cmap="coolwarm",
    s=80,
    edgecolors="black"
)

plt.axhline(-np.log10(0.25), color="gray", linestyle="--", linewidth=1)
plt.title("NES vs -log10(FDR q-val) of GSEA Terms")
plt.xlabel("Normalized Enrichment Score (NES)")
plt.ylabel("-log10(FDR q-val)")
plt.grid(True)

# Label top 5 NES terms
# Label top 5 NES terms (spread out vertically)
top5 = df.nlargest(5, "NES")

for i, row in top5.iterrows():
    plt.text(
        row["NES"] + 0.1,                      # horizontal shift
        -np.log10(row["FDR q-val"] + 1e-8) + (i * 0.2),  # vertical stagger
        row["Term"],
        fontsize=9
    )

plt.tight_layout()
plt.savefig("/Users/bhavanaramachandran/Desktop/gsea_results/hallmark/NES_vs_FDR_dotplot.png")
plt.show()

