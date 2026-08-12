import os, numpy as np, pandas as pd, scanpy as sc
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, roc_auc_score
from scipy.stats import spearmanr
ROOT="/scratch/chahal.d/Spatial_GenAI"; OUT=os.path.join(ROOT,"HE_FoundationModel")
CELL=["Cancer Epithelial","T-cells","Myeloid","CAFs","Endothelial","B-cells","Plasmablasts","PVL","Normal Epithelial"]
TLS=["CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS"]
d=np.load(os.path.join(OUT,"phikon_embeddings.npz"),allow_pickle=True)
X=d["X"]; samp=d["sample"].astype(str); keys=np.array([f"{s}_{b}" for s,b in zip(samp,d["barcode"])])
vm=sc.read_h5ad(os.path.join(ROOT,"cell2location","results","vis_mapped.h5ad"))
comp=vm.obs[CELL].div(vm.obs[CELL].sum(axis=1),axis=0)
ge=sc.read_h5ad(os.path.join(ROOT,"cell2location","data","vis.h5ad")); ge.var_names=[str(x) for x in ge.var_names]
sc.pp.normalize_total(ge,target_sum=1e4); sc.pp.log1p(ge)
present=[g for g in TLS if g in ge.var_names]; sc.tl.score_genes(ge,present,score_name="TLS"); tls=ge.obs["TLS"]
common=set(comp.index)&set(tls.index); m=np.array([k in common for k in keys])
X=X[m]; keys=keys[m]; samp=samp[m]; Y=comp.loc[keys,CELL].to_numpy(); ytls=tls.loc[keys].to_numpy()
Pc=np.zeros_like(Y); Pt=np.zeros_like(ytls)
for s in sorted(set(samp)):
    idx=np.where(samp==s)[0]
    for tr,te in KFold(5,shuffle=True,random_state=0).split(idx):
        tri,tei=idx[tr],idx[te]; sca=StandardScaler().fit(X[tri])
        Pc[tei]=Ridge(alpha=2000).fit(sca.transform(X[tri]),Y[tri]).predict(sca.transform(X[tei]))
        Pt[tei]=Ridge(alpha=2000).fit(sca.transform(X[tri]),ytls[tri]).predict(sca.transform(X[tei]))
rows=[dict(target=c,r2=round(r2_score(Y[:,j],Pc[:,j]),3),spearman=round(spearmanr(Y[:,j],Pc[:,j]).correlation,3)) for j,c in enumerate(CELL)]
tbin=(ytls>=np.median(ytls)).astype(int)
rows.append(dict(target="TLS_score",r2=round(r2_score(ytls,Pt),3),spearman=round(spearmanr(ytls,Pt).correlation,3),auc=round(roc_auc_score(tbin,Pt),3)))
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"he_prediction_metrics_withinslide.csv"),index=False)
print("WITHIN-SLIDE 5-fold CV:\n", res.to_string(index=False)); print("WITHIN_DONE")

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
coords=pd.DataFrame(vm.obsm["spatial"],index=vm.obs_names).loc[keys].to_numpy()
for s in ["CID4465","1160920F"]:
    te=samp==s; xy=coords[te]
    fig,ax=plt.subplots(1,2,figsize=(12,5))
    for a,(vals,t) in zip(ax,[(ytls[te],"TLS — transcriptomics"),(Pt[te],"TLS — predicted from H&E (Phikon)")]):
        p=a.scatter(xy[:,0],-xy[:,1],c=vals,s=9,cmap="magma"); a.set_title(f"{s}  ·  {t}",fontsize=10); a.set_xticks([]); a.set_yticks([]); a.set_aspect("equal"); fig.colorbar(p,ax=a,fraction=0.046)
    fig.tight_layout(); fig.savefig(os.path.join(OUT,f"HE_TLS_within_{s}.png"),dpi=150); plt.close(fig)
print("PLOTS_DONE")
