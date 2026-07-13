# Project spatial-derived signatures onto TCGA-BRCA bulk RNA-seq and test overall
# survival. Signatures: TLS (Cabrita 2020), and PROGENy JAK-STAT & TGFb (Schubert 2018)
# — the axes we found spatially (TLS/lymphoid niche, JAK-STAT in lymphoid niche, TGFb in
# CAF niche). Cohort: UCSC Xena TCGA-BRCA (HiSeqV2, log2 RSEM). See REFERENCES.md.
import os, numpy as np, pandas as pd, anndata as ad, decoupler as dc
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; SV=os.path.join(ROOT,"Survival"); OUT=os.path.join(SV,"results"); os.makedirs(OUT,exist_ok=True)
TLS_GENES=["CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS"]

expr=pd.read_csv(os.path.join(SV,"TCGA_BRCA","HiSeqV2.gz"), sep="\t", index_col=0)   # genes x samples
sv=pd.read_csv(os.path.join(SV,"TCGA_BRCA","BRCA_survival.txt"), sep="\t")[["sample","OS","OS.time"]].dropna()
sv=sv[sv["OS.time"]>0]
samples=[s for s in expr.columns if s in set(sv["sample"])]
expr=expr[samples]; sv=sv.set_index("sample").loc[samples]
print("patients with expression+survival:", len(samples))

# gene z-scores across samples -> TLS mean-z
z=expr.sub(expr.mean(1),axis=0).div(expr.std(1)+1e-9,axis=0)
present=[g for g in TLS_GENES if g in z.index]; print("TLS genes present:",present)
tls=z.loc[present].mean(0)

# PROGENy (JAK-STAT, TGFb) via decoupler on samples x genes
adata=ad.AnnData(expr.T.values.astype("float32"))
adata.obs_names=pd.Index([str(x) for x in expr.columns]); adata.var_names=pd.Index([str(x) for x in expr.index])
net=pd.read_csv(os.path.join(ROOT,"CellComm","progeny_human.csv"))
dc.mt.mlm(data=adata, net=net)
key=[k for k in adata.obsm.keys() if k.startswith("score")][0]
paths=pd.DataFrame(adata.obsm[key], index=adata.obs_names)

sig=pd.DataFrame({"TLS":tls.values,"JAK_STAT":paths["JAK-STAT"].values,"TGFb":paths["TGFb"].values}, index=samples)
df=sig.join(sv); df["OS_months"]=df["OS.time"]/30.44
df.to_csv(os.path.join(OUT,"tcga_signature_scores.csv"))

res=[]
for name in ["TLS","JAK_STAT","TGFb"]:
    d=df[[name,"OS","OS_months"]].dropna().copy()
    d["z"]=(d[name]-d[name].mean())/d[name].std()
    cph=CoxPHFitter().fit(d[["z","OS_months","OS"]], "OS_months","OS")
    hr=float(np.exp(cph.params_["z"])); lo,hi_=np.exp(cph.confidence_intervals_.loc["z"].values); p=float(cph.summary.loc["z","p"])
    grp=d[name]>=d[name].median()
    lr=logrank_test(d.loc[grp,"OS_months"],d.loc[~grp,"OS_months"],d.loc[grp,"OS"],d.loc[~grp,"OS"])
    res.append(dict(signature=name,n=len(d),HR_per_SD=hr,HR_CI_low=lo,HR_CI_high=hi_,cox_p=p,logrank_p=lr.p_value))
    fig,ax=plt.subplots(figsize=(6,5)); km=KaplanMeierFitter()
    km.fit(d.loc[grp,"OS_months"],d.loc[grp,"OS"],label=f"{name} high (n={grp.sum()})"); km.plot_survival_function(ax=ax,ci_show=False)
    km.fit(d.loc[~grp,"OS_months"],d.loc[~grp,"OS"],label=f"{name} low (n={(~grp).sum()})"); km.plot_survival_function(ax=ax,ci_show=False)
    ax.set_title(f"TCGA-BRCA OS by {name}\nHR/SD={hr:.2f} (p={p:.3g}), log-rank p={lr.p_value:.3g}")
    ax.set_xlabel("months"); ax.set_ylabel("OS probability"); ax.set_ylim(0,1); fig.tight_layout()
    fig.savefig(os.path.join(OUT,f"KM_{name}.png"),dpi=150); plt.close(fig)
pd.DataFrame(res).to_csv(os.path.join(OUT,"survival_stats.csv"),index=False)
print(pd.DataFrame(res).round(4).to_string(index=False))
print("SURVIVAL_DONE")
