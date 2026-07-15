# cell2location co-location analysis: NMF of the spot x cell-type abundance matrix to
# find co-occurring cell-type "microenvironments" (Kleshchevnikov 2022, the method-native
# downstream). Reuses vis_mapped.h5ad. See REFERENCES.md.
import os, math, numpy as np, pandas as pd, scanpy as sc
from sklearn.decomposition import NMF
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; RES=os.path.join(ROOT,"cell2location","results")
OUT=os.path.join(RES,"NMF_microenvironments"); os.makedirs(OUT,exist_ok=True)
CELL_TYPES=["Cancer Epithelial","T-cells","Myeloid","CAFs","Endothelial","B-cells","Plasmablasts","PVL","Normal Epithelial"]
N_FACT=6   # < 9 cell types -> groups co-occurring types into microenvironments

vis=sc.read_h5ad(os.path.join(RES,"vis_mapped.h5ad"))
A=np.clip(vis.obs[CELL_TYPES].to_numpy(dtype=float),0,None)
nmf=NMF(n_components=N_FACT, init="nndsvda", random_state=0, max_iter=3000)
W=nmf.fit_transform(A); H=nmf.components_
Hn=H/(H.sum(1,keepdims=True)+1e-9)         # ME x cell-type composition
Wn=W/(W.sum(1,keepdims=True)+1e-9)         # spot x ME (normalised)
me=[f"ME{i+1}" for i in range(N_FACT)]
labels=[f"ME{i+1}: "+" + ".join(np.array(CELL_TYPES)[np.argsort(Hn[i])[::-1][:2]]) for i in range(N_FACT)]

pd.DataFrame(Hn,index=me,columns=CELL_TYPES).to_csv(os.path.join(OUT,"ME_celltype_loadings.csv"))
print("Microenvironment cell-type loadings:\n", pd.DataFrame(Hn,index=me,columns=CELL_TYPES).round(2).to_string())

# loadings heatmap
fig,ax=plt.subplots(figsize=(9,5)); im=ax.imshow(Hn,aspect="auto",cmap="viridis")
ax.set_xticks(range(len(CELL_TYPES))); ax.set_xticklabels(CELL_TYPES,rotation=45,ha="right")
ax.set_yticks(range(N_FACT)); ax.set_yticklabels(me)
for i in range(N_FACT):
  for j in range(len(CELL_TYPES)): ax.text(j,i,f"{Hn[i,j]:.2f}",ha="center",va="center",color="w",fontsize=7)
fig.colorbar(im,fraction=0.03); ax.set_title("cell2location NMF microenvironments — cell-type loadings")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"microenvironment_loadings.png"),dpi=150); plt.close(fig)

# spatial maps of each ME per sample
samples=list(pd.unique(vis.obs["sample"])); samp=vis.obs["sample"].astype(str).to_numpy(); coords=np.asarray(vis.obsm["spatial"],float)
for s in samples:
  m=np.where(samp==s)[0]; ncol=3; nrow=math.ceil(N_FACT/ncol)
  fig,axes=plt.subplots(nrow,ncol,figsize=(4*ncol,4*nrow)); axes=np.atleast_1d(axes).ravel()
  for i in range(N_FACT):
    ax=axes[i]; p=ax.scatter(coords[m,0],-coords[m,1],c=Wn[m,i],s=6,cmap="magma")
    ax.set_title(labels[i],fontsize=8); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal"); fig.colorbar(p,ax=ax,fraction=0.046)
  for j in range(N_FACT,len(axes)): axes[j].axis("off")
  fig.suptitle(f"{s} — NMF microenvironments"); fig.tight_layout(); fig.savefig(os.path.join(OUT,f"ME_spatial_{s}.png"),dpi=140); plt.close(fig)

vis.obs["microenv"]=[me[i] for i in W.argmax(1)]
vis.obs[["sample","microenv"]].to_csv(os.path.join(OUT,"spot_microenvironment_assignment.csv"))
print("dominant-ME counts by sample:\n", pd.crosstab(vis.obs["sample"],vis.obs["microenv"]).to_string())
print("NMF_ME_DONE")
