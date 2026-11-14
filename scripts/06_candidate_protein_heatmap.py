from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

Path("figures").mkdir(exist_ok=True)

matrix = pd.read_csv(
    "data/processed/cll_protein_log2_ratios.csv"
)

results = pd.read_csv(
    "results/candidate_differential_proteins.csv"
)

metadata = pd.read_csv(
    "data/processed/sample_metadata.csv",
    dtype={"sample_id": str}
)

# Exclude the published outlier for primary analysis
samples = metadata.loc[
    ~metadata["published_outlier"],
    "sample_id"
].tolist()

# Select top 15 upregulated and 15 downregulated
up = (
    results[results["direction"] == "Upregulated"]
    .nlargest(15, "mean_log2_ratio")
)

down = (
    results[results["direction"] == "Downregulated"]
    .nsmallest(15, "mean_log2_ratio")
)

selected = pd.concat([up, down])

# Match by UniProt protein-group accession
heatmap_data = (
    matrix.set_index("protein_accession")
    .loc[selected["protein_accession"], samples]
)

# Use gene symbols as row labels
heatmap_data.index = selected["gene_symbol"].values

sns.set_theme(style="white")

plt.figure(figsize=(13, 11))

sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",
    center=0,
    vmin=-3,
    vmax=3,
    linewidths=0.2,
    cbar_kws={
        "label": "log2(CLL / healthy B-cell ratio)"
    }
)

plt.xlabel("CLL patient")
plt.ylabel("Protein")
plt.title(
    "Top 30 Differentially Regulated Proteins in CLL"
)
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()

plt.savefig(
    "figures/04_top30_protein_heatmap.png",
    dpi=300
)
plt.close()

heatmap_data.to_csv(
    "results/top30_protein_abundance_matrix.csv"
)

# Quantify direction consistency
consistency = []

for accession, row in zip(
    selected["protein_accession"],
    heatmap_data.to_numpy()
):
    mean_value = np.mean(row)

    if mean_value > 0:
        fraction = np.mean(row > 0)
    else:
        fraction = np.mean(row < 0)

    consistency.append(fraction)

selected = selected.copy()
selected["direction_consistency"] = consistency

selected.to_csv(
    "results/top30_protein_consistency.csv",
    index=False
)

print("=== HEATMAP RESULTS ===")
print("Proteins:", heatmap_data.shape[0])
print("Patients:", heatmap_data.shape[1])

print("\n=== DIRECTION CONSISTENCY ===")
print(
    selected[[
        "gene_symbol",
        "mean_log2_ratio",
        "direction_consistency"
    ]].to_string(index=False)
)

print("\nSaved:")
print("figures/04_top30_protein_heatmap.png")
print("results/top30_protein_consistency.csv")
