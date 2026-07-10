"""Stage 0: build AnnData objects for cell2location.

Outputs (in cell2location/data/):
  ref.h5ad       reference scRNA-seq, subsampled, raw counts, obs.celltype_major
  vis.h5ad       6 Visium samples concatenated, raw counts, obs.sample, obsm.spatial
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from c2l_io import load_10x, load_spatial_coords

ROOT = "/scratch/chahal.d/Spatial_GenAI"
OUT = os.path.join(ROOT, "cell2location", "data")
os.makedirs(OUT, exist_ok=True)

SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
SUBTYPE = {"1142243F": "TNBC", "1160920F": "TNBC", "CID4290": "ER+",
           "CID4465": "TNBC", "CID44971": "TNBC", "CID4535": "ER+"}

MAX_PER_TYPE = 2500          # cap cells per celltype_major for tractable training
np.random.seed(0)

# ---------------- Reference ----------------
print("Loading reference...", flush=True)
ref_dir = os.path.join(ROOT, "Wu_etal_2021_BRCA_scRNASeq")
ref = load_10x(ref_dir, feature_col=0)

meta = pd.read_csv(os.path.join(ref_dir, "metadata.csv"), index_col=0)
common = ref.obs_names.intersection(meta.index)
print(f"  ref cells={ref.n_obs}, genes={ref.n_vars}; metadata match={len(common)}", flush=True)
ref = ref[common].copy()
meta = meta.loc[common]
for col in ["celltype_major", "orig.ident", "subtype"]:
    if col in meta.columns:
        ref.obs[col] = meta[col].values

# subsample per major cell type
idx = []
for ct, sub in ref.obs.groupby("celltype_major"):
    take = sub.index if len(sub) <= MAX_PER_TYPE else np.random.choice(sub.index, MAX_PER_TYPE, replace=False)
    idx.extend(list(take))
ref = ref[idx].copy()
print("  subsampled reference:", ref.shape, flush=True)
print(ref.obs["celltype_major"].value_counts().to_string(), flush=True)

# basic gene filter
sc.pp.filter_genes(ref, min_cells=5)
ref.write(os.path.join(OUT, "ref.h5ad"))
print("  wrote ref.h5ad", ref.shape, flush=True)

# ---------------- Spatial ----------------
print("Loading spatial samples...", flush=True)
vis_list = []
for s in SAMPLES:
    a = load_10x(os.path.join(ROOT, "samples", s, "filtered_feature_bc_matrix"), feature_col=0)
    coords = load_spatial_coords(os.path.join(ROOT, "samples", s, "spatial"))
    inter = a.obs_names.intersection(coords.index)
    a = a[inter].copy()
    a.obsm["spatial"] = coords.loc[inter, ["pxl_col", "pxl_row"]].values.astype(float)
    a.obs["sample"] = s
    a.obs["subtype"] = SUBTYPE[s]
    a.obs_names = [f"{s}_{b}" for b in a.obs_names]
    print(f"  {s}: {a.shape}", flush=True)
    vis_list.append(a)

vis = ad.concat(vis_list, join="outer", label="batch", index_unique=None)
vis.obs_names_make_unique()
vis.write(os.path.join(OUT, "vis.h5ad"))
print("  wrote vis.h5ad", vis.shape, flush=True)
print("PREP_DONE", flush=True)
