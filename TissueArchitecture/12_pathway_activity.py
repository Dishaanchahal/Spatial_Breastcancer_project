# Pathway activity per spatial niche and per subtype: PROGENy (Schubert 2018)
# via decoupleR (Badia-i-Mompel 2022). See REFERENCES.md.
import os, numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; OUT=os.path.join(ROOT,"TissueArchitecture")
net=pd.read_csv(os.path.join(ROOT,"CellComm","progeny_human.csv"))

ge=sc.read_h5ad(os.path.join(ROOT,"cell2location","data","vis.h5ad"))
# coerce indices to plain object dtype (decoupler 2.x rejects pandas StringDtype)
ge.obs_names=pd.Index([str(x) for x in ge.obs_names]); ge.var_names=pd.Index([str(x) for x in ge.var_names])
sc.pp.normalize_total(ge,target_sum=1e4); sc.pp.log1p(ge)

dc.mt.mlm(data=ge, net=net)                      # decoupler 2.x
key=[k for k in ge.obsm.keys() if k.startswith("score")]
key=key[0] if key else [k for k in ge.obsm.keys() if "mlm" in k and "padj" not in k][0]
obj=ge.obsm[key]
scores=obj.copy() if isinstance(obj,pd.DataFrame) else pd.DataFrame(np.asarray(obj), index=ge.obs_names)
scores.index=ge.obs_names
scores.to_csv(os.path.join(OUT,"progeny_scores_per_spot.csv"))
print("pathways:", list(scores.columns))

nl=pd.read_csv(os.path.join(OUT,"niche_labels.csv")); nl["barcode"]=nl["barcode"].astype(str)
nl=nl.set_index("barcode").reindex(ge.obs_names)
byniche=scores.groupby(nl["niche"].values).mean(); byniche.to_csv(os.path.join(OUT,"progeny_by_niche.csv"))
bysub=scores.groupby(ge.obs["subtype"].astype(str).values).mean(); bysub.to_csv(os.path.join(OUT,"progeny_by_subtype.csv"))

def heat(df,title,fn):
    fig,ax=plt.subplots(figsize=(max(8,0.55*df.shape[1]),0.6*df.shape[0]+2))
    v=float(np.nanmax(np.abs(df.values))) or 1.0
    im=ax.imshow(df.values,aspect="auto",cmap="RdBu_r",vmin=-v,vmax=v)
    ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns,rotation=90,fontsize=8)
    ax.set_yticks(range(df.shape[0])); ax.set_yticklabels([str(i) for i in df.index])
    ax.set_title(title); fig.colorbar(im,label="PROGENy activity (MLM)"); fig.tight_layout()
    fig.savefig(os.path.join(OUT,fn),dpi=150); plt.close(fig)
byniche.index=[f"niche {i}" for i in byniche.index]
heat(byniche,"PROGENy pathway activity by spatial niche","progeny_by_niche_heatmap.png")
heat(bysub,"PROGENy pathway activity by subtype","progeny_by_subtype_heatmap.png")
print("PATHWAY_DONE")
