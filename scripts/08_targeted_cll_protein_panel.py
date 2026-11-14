from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

Path("figures").mkdir(exist_ok=True)

df = pd.read_csv(
    "results/cll_bcell_targeted_protein_panel.csv"
)

df = df.sort_values("mean_log2_ratio")

colors = np.where(
    df["high_priority"],
    "#477FBD",
    "#B0B0B0"
)

fig, ax = plt.subplots(figsize=(10, 8))

bars = ax.barh(
    df["gene_symbol"],
    df["mean_log2_ratio"],
    color=colors
)

ax.axvline(0, color="black", linewidth=1)

ax.axvline(
    -0.5, color="grey", linestyle="--", alpha=0.6
)
ax.axvline(
    0.5, color="grey", linestyle="--", alpha=0.6
)

for bar, (_, row) in zip(bars, df.iterrows()):
    value = row["mean_log2_ratio"]
    offset = 0.06 if value >= 0 else -0.06

    ax.text(
        value + offset,
        bar.get_y() + bar.get_height() / 2,
        f"{value:+.2f}",
        va="center",
        ha="left" if value >= 0 else "right",
        fontsize=9
    )

ax.set_xlim(-1.8, 3.0)
ax.set_xlabel(
    "Mean log2(CLL / healthy B-cell ratio)"
)
ax.set_ylabel("Protein")
ax.set_title(
    "CLL-Associated B-Cell and Signalling Proteins"
)

from matplotlib.patches import Patch

ax.legend(
    handles=[
        Patch(
            facecolor="#477FBD",
            label="Meets exploratory priority criteria"
        ),
        Patch(
            facecolor="#B0B0B0",
            label="Does not meet all criteria"
        )
    ],
    loc="lower right"
)

plt.tight_layout()

plt.savefig(
    "figures/05_targeted_cll_protein_panel.png",
    dpi=300
)
plt.close()

print("Targeted CLL protein panel saved.")
print("Proteins plotted:", len(df))
print(
    "High-priority proteins:",
    df["high_priority"].sum()
)
