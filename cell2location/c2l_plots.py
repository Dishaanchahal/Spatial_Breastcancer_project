"""Stage 3: figures + composition table from cell2location posteriors.

Outputs (in cell2location/results/):
  abundance_<SAMPLE>.png       spatial maps of every cell type (q05 abundance)
  tme_composition_c2l.png      per-sample composition, faceted by subtype
  tme_composition_c2l.csv      mean-abundance-derived proportions
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/scratch/chahal.d/Spatial_GenAI"
RES = os.path.join(ROOT, "cell2location", "results")

vis = sc.read_h5ad(os.path.join(RES, "vis_mapped.h5ad"))
cell_types = list(vis.uns["mod"]["factor_names"])
samples = list(pd.unique(vis.obs["sample"]))
subtype = dict(zip(vis.obs["sample"], vis.obs["subtype"]))

# ---- per-sample spatial abundance maps ----
for s in samples:
    sub = vis[vis.obs["sample"] == s]
    xy = sub.obsm["spatial"]
    n = len(cell_types)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for i, ct in enumerate(cell_types):
        ax = axes[i]
        vals = np.asarray(sub.obs[ct]).astype(float)
        sc_ = ax.scatter(xy[:, 0], -xy[:, 1], c=vals, s=6, cmap="magma")
        ax.set_title(ct, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        fig.colorbar(sc_, ax=ax, fraction=0.046)
    for j in range(n, len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"{s} ({subtype[s]}) — cell2location q05 abundance", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, f"abundance_{s}.png"), dpi=150)
    plt.close(fig)
    print("wrote", f"abundance_{s}.png", flush=True)

# ---- composition table (mean abundance -> proportion) ----
df = vis.obs[["sample"] + cell_types].groupby("sample").mean()
prop = df.div(df.sum(axis=1), axis=0)
prop["subtype"] = [subtype[s] for s in prop.index]
prop.to_csv(os.path.join(RES, "tme_composition_c2l.csv"))
print("composition (proportions):", flush=True)
print(prop.round(3).to_string(), flush=True)

# ---- stacked bar, faceted by subtype ----
order = sorted(prop.index, key=lambda s: (prop.loc[s, "subtype"], s))
fig, ax = plt.subplots(figsize=(10, 6))
bottom = np.zeros(len(order))
cmap = plt.get_cmap("tab10")
for k, ct in enumerate(cell_types):
    vals = prop.loc[order, ct].values.astype(float)
    ax.bar(range(len(order)), vals, bottom=bottom, label=ct, color=cmap(k % 10))
    bottom += vals
ax.set_xticks(range(len(order)))
ax.set_xticklabels([f"{s}\n{prop.loc[s,'subtype']}" for s in order], fontsize=8)
ax.set_ylabel("Proportion (cell2location)")
ax.set_title("TME composition — cell2location deconvolution")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(RES, "tme_composition_c2l.png"), dpi=150)
plt.close(fig)
print("wrote tme_composition_c2l.png", flush=True)
print("PLOTS_DONE", flush=True)
