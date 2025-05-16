import gseapy as gp
import pandas as pd
import os

# === Step 1: Load ranked gene list ===
rnk_path = "/Users/bhavanaramachandran/Desktop/gene_list_for_gsea.csv"
rnk = pd.read_csv(rnk_path)

# Confirm required columns exist
assert {'gene_name', 'stat'}.issubset(rnk.columns), "CSV must have 'gene_name' and 'stat' columns"

# gseapy requires a 2-column DataFrame: gene name and ranking score
rnk = rnk[['gene_name', 'stat']]

# === Step 2: Define gene sets ===
gene_sets = {
    "hallmark": "/Users/bhavanaramachandran/Desktop/h.all.v2024.1.Hs.symbols.gmt",
    "canonical": "/Users/bhavanaramachandran/Desktop/c2.cp.v2024.1.Hs.symbols.gmt"
}

# === Step 3: Output directory ===
outdir_base = "/Users/bhavanaramachandran/Desktop/gsea_results"
os.makedirs(outdir_base, exist_ok=True)

# === Step 4: Run GSEA for each gene set ===
for name, gmt_path in gene_sets.items():
    print(f"Running GSEA for: {name}")
    gp.prerank(
        rnk=rnk,
        gene_sets=gmt_path,
        outdir=os.path.join(outdir_base, name),
        min_size=15,
        max_size=500,
        permutation_num=1000,  # Increase for publication-level results
        seed=42,
        verbose=True
    )

