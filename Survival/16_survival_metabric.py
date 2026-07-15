# METABRIC survival validation (independent, well-powered cohort; Curtis 2012 / Pereira 2016,
# obtained via a GitHub mirror of the cBioPortal study). Projects TLS + PROGENy JAK-STAT/TGFb
# signatures and tests OS overall + within Basal / ER-HER2- subsets. See REFERENCES.md.
import os, re, numpy as np, pandas as pd, anndata as ad, decoupler as dc
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; SV=os.path.join(ROOT,"Survival"); MB=os.path.join(SV,"METABRIC"); OUT=os.path.join(SV,"results")
TLS_GENES=["CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS"]
norm=lambda s: re.sub(r"[^A-Za-z0-9]","",str(s)).upper()

expr=pd.read_csv(os.path.join(MB,"metabric_expr.tsv.gz"), sep="\t", index_col=0)  # genes x samples
clin=pd.read_csv(os.path.join(MB,"metabric_clin.csv"))
clin["key"]=clin["PATIENT_ID"].map(norm); clin=clin.drop_duplicates("key").set_index("key")
expr.columns=[norm(c) for c in expr.columns]
samples=[c for c in expr.columns if c in clin.index]
expr=expr[samples]; clin=clin.loc[samples]
print("METABRIC samples with expr+clin:", len(samples))

# signatures
z=expr.sub(expr.mean(1),axis=0).div(expr.std(1)+1e-9,axis=0)
present=[g for g in TLS_GENES if g in z.index]; print("TLS genes present:",present)
tls=z.loc[present].mean(0)
adata=ad.AnnData(expr.T.values.astype("float32")); adata.obs_names=pd.Index(expr.columns); adata.var_names=pd.Index([str(x) for x in expr.index])
net=pd.read_csv(os.path.join(ROOT,"CellComm","progeny_human.csv"))
dc.mt.mlm(data=adata, net=net); key=[k for k in adata.obsm if k.startswith("score")][0]
paths=pd.DataFrame(adata.obsm[key], index=adata.obs_names)

df=pd.DataFrame({"TLS":tls,"JAK_STAT":paths["JAK-STAT"],"TGFb":paths["TGFb"]})
df["OS_MONTHS"]=pd.to_numeric(clin["OS_MONTHS"],errors="coerce")
df["event"]=clin["OS_STATUS"].astype(str).str.contains("DECEASED",case=False).astype(int)
df["CLAUDIN"]=clin["CLAUDIN_SUBTYPE"]; df["THREEGENE"]=clin["THREEGENE"]
df=df[df["OS_MONTHS"]>0]
df.to_csv(os.path.join(OUT,"metabric_signature_scores.csv"))

groups={"ALL":df.index,
        "Basal":df.index[df["CLAUDIN"]=="Basal"],
        "ER-HER2-":df.index[df["THREEGENE"].astype(str).str.contains("ER-/HER2-",na=False)]}
res=[]
for g,idx in groups.items():
    b=df.loc[idx]
    for name in ["TLS","JAK_STAT","TGFb"]:
        d=b[[name,"OS_MONTHS","event"]].dropna().copy(); ev=int(d["event"].sum())
        if len(d)<20 or ev<5: res.append(dict(group=g,signature=name,n=len(d),events=ev,note="underpowered")); continue
        d["z"]=(d[name]-d[name].mean())/d[name].std()
        cph=CoxPHFitter().fit(d[["z","OS_MONTHS","event"]],"OS_MONTHS","event")
        hr=float(np.exp(cph.params_["z"])); lo,hi=np.exp(cph.confidence_intervals_.loc["z"].values); p=float(cph.summary.loc["z","p"])
        grp=d[name]>=d[name].median()
        lr=logrank_test(d.loc[grp,"OS_MONTHS"],d.loc[~grp,"OS_MONTHS"],d.loc[grp,"event"],d.loc[~grp,"event"])
        res.append(dict(group=g,signature=name,n=len(d),events=ev,HR_per_SD=round(hr,3),CI_low=round(lo,3),CI_high=round(hi,3),cox_p=round(p,5),logrank_p=round(lr.p_value,5)))
        if g in ("ALL","ER-HER2-"):
            fig,ax=plt.subplots(figsize=(6,5)); km=KaplanMeierFitter()
            km.fit(d.loc[grp,"OS_MONTHS"],d.loc[grp,"event"],label=f"{name} high (n={int(grp.sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
            km.fit(d.loc[~grp,"OS_MONTHS"],d.loc[~grp,"event"],label=f"{name} low (n={int((~grp).sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
            ax.set_title(f"METABRIC {g} — OS by {name}\nHR/SD={hr:.2f} p={p:.2g}, log-rank {lr.p_value:.2g}"); ax.set_ylim(0,1); ax.set_xlabel("months")
            fig.tight_layout(); fig.savefig(os.path.join(OUT,f"METABRIC_KM_{g}_{name}.png"),dpi=150); plt.close(fig)
out=pd.DataFrame(res); out.to_csv(os.path.join(OUT,"survival_metabric_stats.csv"),index=False)
print(out.to_string(index=False)); print("METABRIC_SURVIVAL_DONE")
