"""Stage 1+2: cell2location reference signatures + spatial mapping (GPU).

Follows the official cell2location workflow:
  1) RegressionModel on the scRNA-seq reference -> per-cell-type signatures.
  2) Cell2location model maps those signatures onto the Visium spots.

Outputs (in cell2location/results/):
  inf_aver.csv                per-gene, per-cell-type reference signature
  vis_mapped.h5ad             Visium AnnData with posterior cell abundances
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import cell2location
from cell2location.models import RegressionModel, Cell2location

ROOT = "/scratch/chahal.d/Spatial_GenAI"
DATA = os.path.join(ROOT, "cell2location", "data")
RES = os.path.join(ROOT, "cell2location", "results")
os.makedirs(RES, exist_ok=True)

print("torch", torch.__version__, "cuda_available", torch.cuda.is_available(), flush=True)
if torch.cuda.is_available():
    print("gpu:", torch.cuda.get_device_name(0), flush=True)

# ------------------------------------------------------------------ reference
ref = sc.read_h5ad(os.path.join(DATA, "ref.h5ad"))
batch_key = "orig.ident" if "orig.ident" in ref.obs.columns else None
RegressionModel.setup_anndata(adata=ref, batch_key=batch_key, labels_key="celltype_major")
reg = RegressionModel(ref)
reg.view_anndata_setup()
reg.train(max_epochs=250, accelerator="gpu" if torch.cuda.is_available() else "cpu")
ref = reg.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})

# per-cell-type expression signature
if "means_per_cluster_mu_fg" in ref.varm.keys():
    inf_aver = ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{c}" for c in ref.uns["mod"]["factor_names"]]].copy()
else:
    inf_aver = ref.var[
        [f"means_per_cluster_mu_fg_{c}" for c in ref.uns["mod"]["factor_names"]]].copy()
inf_aver.columns = ref.uns["mod"]["factor_names"]
inf_aver.to_csv(os.path.join(RES, "inf_aver.csv"))
print("reference signatures:", inf_aver.shape, flush=True)

# ------------------------------------------------------------------ spatial
vis = sc.read_h5ad(os.path.join(DATA, "vis.h5ad"))
shared = [g for g in vis.var_names if g in inf_aver.index]
print(f"shared genes vis∩ref = {len(shared)}", flush=True)
vis = vis[:, shared].copy()
inf_aver = inf_aver.loc[shared, :]

Cell2location.setup_anndata(adata=vis, batch_key="sample")
mod = Cell2location(
    vis, cell_state_df=inf_aver,
    N_cells_per_location=30,   # prior mean cells/spot for tumor Visium
    detection_alpha=20,        # recommended for within-slide normalisation
)
mod.view_anndata_setup()
mod.train(max_epochs=20000,
          accelerator="gpu" if torch.cuda.is_available() else "cpu",
          batch_size=None, train_size=1)

vis = mod.export_posterior(vis, sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs})

# put posterior 5% quantile abundances into obs (conservative estimate)
vis.obs[vis.uns["mod"]["factor_names"]] = vis.obsm["q05_cell_abundance_w_sf"]
vis.write(os.path.join(RES, "vis_mapped.h5ad"))
print("wrote vis_mapped.h5ad", vis.shape, flush=True)
print("C2L_RUN_DONE", flush=True)
