from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

Path("figures").mkdir(exist_ok=True)

df = pd.read_csv("data/processed/cll_protein_log2_ratios.csv")
metadata = pd.read_csv(
    "data/processed/sample_metadata.csv",
    dtype={"sample_id": str}
)

samples = metadata["sample_id"].tolist()
X = df[samples].T

sns.set_theme(style="whitegrid")

# Figure 1: Abundance distributions
plt.figure(figsize=(12, 6))

sns.boxplot(
    data=X.T,
    color="skyblue",
    showfliers=False
)

plt.axhline(0, color="red", linestyle="--", linewidth=1)
plt.xlabel("CLL patient")
plt.ylabel("Protein abundance (log2 ratio)")
plt.title("CLL Proteomics: Sample Abundance Distributions")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(
    "figures/01_protein_abundance_boxplot.png",
    dpi=300
)
plt.close()

# Figure 2: PCA
# Centre and scale each protein across patients
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame({
    "sample_id": samples,
    "PC1": pcs[:, 0],
    "PC2": pcs[:, 1],
    "batch": metadata["tmt_batch"].values,
    "outlier": metadata["published_outlier"].values
})

plt.figure(figsize=(9, 7))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="batch",
    style="outlier",
    s=130
)

for _, row in pca_df.iterrows():
    plt.annotate(
        row["sample_id"],
        (row["PC1"], row["PC2"]),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=9
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
)
plt.title("PCA of CLL Proteomic Profiles")
plt.tight_layout()

plt.savefig(
    "figures/02_proteomics_pca.png",
    dpi=300
)
plt.close()

pca_df.to_csv(
    "data/processed/pca_coordinates.csv",
    index=False
)

print("=== QC RESULTS ===")
print("Proteins analysed:", X.shape[1])
print("Patients analysed:", X.shape[0])

print("\nPCA explained variance:")
print(
    f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}%"
)
print(
    f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}%"
)

print("\nPCA coordinates:")
print(pca_df.round(2).to_string(index=False))

print("\nFigures saved:")
print("figures/01_protein_abundance_boxplot.png")
print("figures/02_proteomics_pca.png")
