# Cross-validation: map the SAME reference signatures (inf_aver.csv) onto an
# INDEPENDENT 10x Visium breast-cancer dataset (Block A, Sections 1-2).
import os, sys, numpy as np, pandas as pd, scanpy as sc, anndata as ad, torch
import cell2location
from cell2location.models import Cell2location
ROOT="/scratch/chahal.d/Spatial_GenAI"; XV=os.path.join(ROOT,"ExternalValidation")
OUT=os.path.join(XV,"results"); os.makedirs(OUT,exist_ok=True)
sys.path.insert(0, os.path.join(ROOT,"cell2location"))
from c2l_io import load_10x, load_spatial_coords
SAMPLES=["V1_Breast_Cancer_Block_A_Section_1","V1_Breast_Cancer_Block_A_Section_2"]

vis_list=[]
for s in SAMPLES:
    a=load_10x(os.path.join(XV,"samples",s,"filtered_feature_bc_matrix"), feature_col=1)  # 10x: gene_symbol col
    coords=load_spatial_coords(os.path.join(XV,"samples",s,"spatial"))
    inter=a.obs_names.intersection(coords.index); a=a[inter].copy()
    a.obsm["spatial"]=coords.loc[inter,["pxl_col","pxl_row"]].values.astype(float)
    a.obs["sample"]=s; a.obs_names=[f"{s}_{b}" for b in a.obs_names]
    print(s, a.shape, flush=True); vis_list.append(a)
vis=ad.concat(vis_list, join="outer", label="batch", index_unique=None); vis.obs_names_make_unique()
vis.write(os.path.join(OUT,"vis_xval.h5ad"))                 # full genes for TLS/PROGENy

inf=pd.read_csv(os.path.join(ROOT,"cell2location","results","inf_aver.csv"), index_col=0)
shared=[g for g in vis.var_names if g in inf.index]
print("shared genes vis∩reference:", len(shared), flush=True)
vm=vis[:, shared].copy(); inf=inf.loc[shared]
Cell2location.setup_anndata(adata=vm, batch_key="sample")
mod=Cell2location(vm, cell_state_df=inf, N_cells_per_location=30, detection_alpha=20)
mod.train(max_epochs=20000, accelerator="gpu" if torch.cuda.is_available() else "cpu", batch_size=None, train_size=1)
vm=mod.export_posterior(vm, sample_kwargs={"num_samples":1000,"batch_size":vm.n_obs})
vm.obs[vm.uns["mod"]["factor_names"]]=vm.obsm["q05_cell_abundance_w_sf"]
vm.write(os.path.join(OUT,"vis_xval_mapped.h5ad"))
print("XVAL_C2L_DONE", flush=True)
