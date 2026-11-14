from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

Path("figures").mkdir(exist_ok=True)

df = pd.read_csv(
    "results/all_protein_regulation_results.csv"
)

df["status"] = "Not selected"

df.loc[
    (df["fdr"] < 0.05) &
    (df["mean_log2_ratio"] >= 0.5),
    "status"
] = "Upregulated"

df.loc[
    (df["fdr"] < 0.05) &
    (df["mean_log2_ratio"] <= -0.5),
    "status"
] = "Downregulated"

sns.set_theme(style="whitegrid")

plt.figure(figsize=(11, 8))

palette = {
    "Not selected": "#B0B0B0",
    "Upregulated": "#D97757",
    "Downregulated": "#477FBD"
}

for status in [
    "Not selected",
    "Downregulated",
    "Upregulated"
]:
    subset = df[df["status"] == status]

    plt.scatter(
        subset["mean_log2_ratio"],
        subset["neg_log10_fdr"],
        c=palette[status],
        s=13,
        alpha=0.65,
        label=f"{status} (n={len(subset)})",
        rasterized=True
    )

plt.axvline(-0.5, color="grey", linestyle="--")
plt.axvline(0.5, color="grey", linestyle="--")
plt.axhline(
    -np.log10(0.05),
    color="grey",
    linestyle="--"
)

# Label strongest candidates by absolute effect size
top = df[df["candidate"]].nlargest(
    12, "abs_mean_log2_ratio"
)

for _, row in top.iterrows():
    plt.annotate(
        row["gene_symbol"],
        (
            row["mean_log2_ratio"],
            row["neg_log10_fdr"]
        ),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=8
    )

plt.xlabel("Mean log2(CLL / healthy B-cell ratio)")
plt.ylabel("-log10(published FDR-adjusted p-value)")
plt.title("CLL Proteomics: Differential Protein Regulation")
plt.legend(frameon=True)
plt.tight_layout()

plt.savefig(
    "figures/03_cll_proteomics_volcano.png",
    dpi=300
)
plt.close()

print("Volcano plot saved successfully.")
print("Total proteins:", len(df))
print("Candidate proteins:", df["candidate"].sum())
print("Upregulated:",
      (df["status"] == "Upregulated").sum())
print("Downregulated:",
      (df["status"] == "Downregulated").sum())
