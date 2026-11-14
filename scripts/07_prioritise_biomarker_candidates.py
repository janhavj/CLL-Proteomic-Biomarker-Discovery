from pathlib import Path
import numpy as np
import pandas as pd

Path("results").mkdir(exist_ok=True)

matrix = pd.read_csv(
    "data/processed/cll_protein_log2_ratios.csv"
)

results = pd.read_csv(
    "results/all_protein_regulation_results.csv"
)

metadata = pd.read_csv(
    "data/processed/sample_metadata.csv",
    dtype={"sample_id": str}
)

samples = metadata.loc[
    ~metadata["published_outlier"],
    "sample_id"
].tolist()

# Calculate consistency for every protein
values = matrix[samples].to_numpy()

means = values.mean(axis=1)

consistency = np.where(
    means >= 0,
    (values > 0).mean(axis=1),
    (values < 0).mean(axis=1)
)

quality = matrix[[
    "protein_accession",
    "unique_peptides",
    "sequence_coverage",
    "psms"
]].copy()

quality["direction_consistency"] = consistency

ranked = results.merge(
    quality,
    on="protein_accession",
    how="left",
    validate="one_to_one"
)

# Exploratory prioritisation criteria
ranked["high_priority"] = (
    (ranked["fdr"] < 0.05) &
    (ranked["abs_mean_log2_ratio"] >= 0.5) &
    (ranked["direction_consistency"] >= 0.80) &
    (ranked["unique_peptides"] >= 2)
)

priority = ranked[
    ranked["high_priority"]
].copy()

priority = priority.sort_values(
    ["abs_mean_log2_ratio", "fdr"],
    ascending=[False, True]
)

ranked.to_csv(
    "results/all_protein_prioritisation.csv",
    index=False
)

priority.to_csv(
    "results/high_priority_protein_candidates.csv",
    index=False
)

print("=== CANDIDATE PRIORITISATION ===")
print("Total protein groups:", len(ranked))
print("Initial candidates:", ranked["candidate"].sum())
print("High-priority candidates:", len(priority))

print("\n=== TOP 20 PRIORITISED PROTEINS ===")

columns = [
    "gene_symbol",
    "mean_log2_ratio",
    "fdr",
    "direction_consistency",
    "unique_peptides"
]

print(
    priority[columns]
    .head(20)
    .to_string(index=False)
)

print("\nSaved:")
print("results/all_protein_prioritisation.csv")
print("results/high_priority_protein_candidates.csv")
