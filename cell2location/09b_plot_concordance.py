import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, os
ROOT = "/scratch/chahal.d/Spatial_GenAI"; OUT = os.path.join(ROOT, "Validation_RCTD")
m = pd.read_csv(os.path.join(OUT, "c2l_vs_rctd_merged.csv"))
r = np.corrcoef(m["c2l"], m["rctd"])[0, 1]
cts = list(pd.unique(m["cell_type"])); cmap = plt.get_cmap("tab10")
fig, ax = plt.subplots(figsize=(8, 7))
for i, ct in enumerate(cts):
    d = m[m["cell_type"] == ct]
    ax.scatter(d["c2l"], d["rctd"], s=60, color=cmap(i % 10), label=ct)
lim = max(m["c2l"].max(), m["rctd"].max()) * 1.05
ax.plot([0, lim], [0, lim], "--", color="grey")
ax.set_xlabel("cell2location proportion"); ax.set_ylabel("RCTD proportion")
ax.set_title(f"cell2location vs RCTD composition (Pearson r = {r:.2f})")
ax.legend(fontsize=8, loc="upper left", frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "c2l_vs_rctd_scatter.png"), dpi=150)
print(f"Pearson r = {r:.3f}; wrote scatter")
