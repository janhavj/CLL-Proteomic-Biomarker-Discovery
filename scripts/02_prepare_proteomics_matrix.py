from pathlib import Path
import pandas as pd
import numpy as np

input_file = Path(
    "data/raw/134638_0_supp_38937_p0y7zb.xlsx"
)
output_dir = Path("data/processed")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_excel(input_file, sheet_name="CLL proteome")

# Standardise Excel column headers to strings
df.columns = df.columns.map(str)

sample_columns = [
    "4621*", "3999", "2497", "568", "653",
    "3405", "3156", "2886", "2483", "625",
    "4140", "1194", "1125", "2459"
]

annotation_columns = [
    "Protein group accession",
    "Description",
    "Symbol",
    "ΣCoverage",
    "Σ# Unique Peptides",
    "Σ# PSMs"
]

clean = df[annotation_columns + sample_columns].copy()

clean = clean.rename(columns={
    "Protein group accession": "protein_accession",
    "Description": "description",
    "Symbol": "gene_symbol",
    "ΣCoverage": "sequence_coverage",
    "Σ# Unique Peptides": "unique_peptides",
    "Σ# PSMs": "psms",
    "4621*": "4621"
})

patient_columns = [
    "4621", "3999", "2497", "568", "653",
    "3405", "3156", "2886", "2483", "625",
    "4140", "1194", "1125", "2459"
]

for col in patient_columns:
    clean[col] = pd.to_numeric(clean[col], errors="coerce")

clean.to_csv(
    output_dir / "cll_protein_log2_ratios.csv",
    index=False
)

metadata = pd.DataFrame({
    "sample_id": patient_columns,
    "tmt_batch": ["A"] * 7 + ["B"] * 7,
    "published_outlier": [
        True, False, False, False, False, False, False,
        False, False, False, False, False, False, False
    ]
})

metadata.to_csv(
    output_dir / "sample_metadata.csv",
    index=False
)

print("=== PROTEOMICS MATRIX ===")
print("Proteins:", len(clean))
print("Patient samples:", len(patient_columns))
print("Matrix shape:", clean[patient_columns].shape)

print("\nMissing abundance values:")
print(clean[patient_columns].isna().sum().sum())

print("\n=== SAMPLE METADATA ===")
print(metadata.to_string(index=False))

print("\n=== ABUNDANCE SUMMARY ===")
print(
    clean[patient_columns]
    .describe()
    .loc[["mean", "std", "min", "max"]]
    .round(3)
    .to_string()
)

print("\nFiles saved to:", output_dir)
