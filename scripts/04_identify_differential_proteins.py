from pathlib import Path
import numpy as np
import pandas as pd

Path("results").mkdir(exist_ok=True)

path = "data/raw/134638_0_supp_38937_p0y7zb.xlsx"
df = pd.read_excel(path, sheet_name="CLL proteome")
df.columns = df.columns.map(str)

# The original study excluded patient 4621
samples = [
    "3999", "2497", "568", "653", "3405", "3156",
    "2886", "2483", "625", "4140", "1194",
    "1125", "2459"
]

results = df[[
    "Protein group accession",
    "Symbol",
    "Regulation score",
    "-Log10(p value)"
]].copy()

results.columns = [
    "protein_accession",
    "gene_symbol",
    "regulation_score",
    "neg_log10_fdr"
]

# Recover the published FDR-adjusted p-values
results["fdr"] = 10 ** (-results["neg_log10_fdr"])

# Descriptive mean log2 ratio across 13 patients
results["mean_log2_ratio"] = df[samples].mean(axis=1)

# Magnitude of the average protein change
results["abs_mean_log2_ratio"] = (
    results["mean_log2_ratio"].abs()
)

# Exploratory candidate threshold:
# FDR < 0.05 and |mean log2 ratio| >= 0.5
results["candidate"] = (
    (results["fdr"] < 0.05) &
    (results["abs_mean_log2_ratio"] >= 0.5)
)

results["direction"] = np.where(
    results["mean_log2_ratio"] > 0,
    "Upregulated",
    "Downregulated"
)

results = results.sort_values(
    ["candidate", "abs_mean_log2_ratio"],
    ascending=[False, False]
)

results.to_csv(
    "results/all_protein_regulation_results.csv",
    index=False
)

candidates = results[results["candidate"]].copy()

candidates.to_csv(
    "results/candidate_differential_proteins.csv",
    index=False
)

print("=== DIFFERENTIAL PROTEIN ANALYSIS ===")
print("Total proteins:", len(results))
print("FDR < 0.05:", (results["fdr"] < 0.05).sum())
print("Candidate proteins:", len(candidates))

print("\nUpregulated candidates:",
      (candidates["direction"] == "Upregulated").sum())

print("Downregulated candidates:",
      (candidates["direction"] == "Downregulated").sum())

print("\n=== TOP 15 CANDIDATES ===")
print(
    candidates[[
        "gene_symbol",
        "mean_log2_ratio",
        "regulation_score",
        "fdr",
        "direction"
    ]].head(15).to_string(index=False)
)

print("\nResults saved to results/")
