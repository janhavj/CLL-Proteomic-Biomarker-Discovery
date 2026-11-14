from pathlib import Path
import pandas as pd

data_dir = Path("data/raw")

files = {
    "Table S1 - Clinical metadata": "134638_0_supp_38936_p0ykzb.xlsx",
    "Table S2 - Quantitative proteomics": "134638_0_supp_38937_p0y7zb.xlsx",
    "Table S3 - Molecular subgroups": "134638_0_supp_38938_p0ykzb.xlsx",
}

for label, filename in files.items():
    print("\n" + "=" * 80)
    print(label)
    print("=" * 80)

    path = data_dir / filename
    workbook = pd.ExcelFile(path)

    print("Worksheets:", workbook.sheet_names)

    for sheet in workbook.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet, header=None)

        print(f"\nSheet: {sheet}")
        print("Dimensions:", df.shape)
        print("\nFirst 8 rows and 12 columns:")
        print(df.iloc[:8, :12].to_string(
            index=False, header=False, max_colwidth=24
        ))
