# Cross-validation downstream: composition vs our cohort, TLS signature, PROGENy.
import os, numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; XV=os.path.join(ROOT,"ExternalValidation"); OUT=os.path.join(XV,"results")
CELL_TYPES=["Cancer Epithelial","T-cells","Myeloid","CAFs","Endothelial","B-cells","Plasmablasts","PVL","Normal Epithelial"]
TLS_GENES=["CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS"]

vm=sc.read_h5ad(os.path.join(OUT,"vis_xval_mapped.h5ad"))
samples=vm.obs["sample"].astype(str).values; coords=np.asarray(vm.obsm["spatial"],float)

# --- composition, compared to our cohort mean ---
comp=vm.obs[CELL_TYPES].copy(); comp["sample"]=samples
prop=comp.groupby("sample").mean(); prop=prop.div(prop.sum(1),axis=0)
ours=pd.read_csv(os.path.join(ROOT,"cell2location","results","tme_composition_c2l.csv"))
ours_prop=ours[CELL_TYPES].div(ours[CELL_TYPES].sum(1),axis=0).mean()
prop.loc["OUR_COHORT_mean"]=ours_prop.values
prop.to_csv(os.path.join(OUT,"xval_composition_vs_ours.csv"))
print("composition (proportions):\n", prop.round(3).to_string())

# --- TLS + PROGENy on full-gene expression ---
ge=sc.read_h5ad(os.path.join(OUT,"vis_xval.h5ad"))
ge.obs_names=pd.Index([str(x) for x in ge.obs_names]); ge.var_names=pd.Index([str(x) for x in ge.var_names])
sc.pp.normalize_total(ge,target_sum=1e4); sc.pp.log1p(ge)
present=[g for g in TLS_GENES if g in ge.var_names]; sc.tl.score_genes(ge,present,score_name="TLS")
print("TLS genes present:",present)
net=pd.read_csv(os.path.join(ROOT,"CellComm","progeny_human.csv"))
dc.mt.mlm(data=ge, net=net)
key=[k for k in ge.obsm.keys() if k.startswith("score")][0]
paths=pd.DataFrame(ge.obsm[key]); paths.index=ge.obs_names
ge=ge[vm.obs_names].copy(); tls=ge.obs["TLS"].to_numpy()
paths=paths.reindex(vm.obs_names)

# --- figures per sample: TLS, JAK-STAT, TGFb ---
def spat(vals,title,fn,cmap="viridis"):
    fig,axes=plt.subplots(1,len(pd.unique(samples)),figsize=(6*len(pd.unique(samples)),5))
    axes=np.atleast_1d(axes)
    for ax,s in zip(axes,pd.unique(samples)):
        m=np.where(samples==s)[0]
        sc_=ax.scatter(coords[m,0],-coords[m,1],c=np.asarray(vals)[m],s=8,cmap=cmap)
        ax.set_title(f"{s}\n{title}",fontsize=9); ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
        fig.colorbar(sc_,ax=ax,fraction=0.046)
    fig.tight_layout(); fig.savefig(os.path.join(OUT,fn),dpi=140); plt.close(fig)
spat(tls,"TLS score (Cabrita)","xval_TLS.png")
spat(paths["JAK-STAT"].to_numpy(),"JAK-STAT (PROGENy)","xval_JAK-STAT.png",cmap="RdBu_r")
spat(paths["TGFb"].to_numpy(),"TGFb (PROGENy)","xval_TGFb.png",cmap="RdBu_r")

# pathway means per sample
pm=paths.copy(); pm["sample"]=samples; pm.groupby("sample").mean().to_csv(os.path.join(OUT,"xval_progeny_by_sample.csv"))
print("XVAL_DOWNSTREAM_DONE")
