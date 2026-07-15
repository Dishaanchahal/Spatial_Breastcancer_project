# Recurrence-endpoint survival on TCGA-BRCA: PFI (progression-free interval) and DFI
# (disease-free interval) — more events than OS. Full cohort + TNBC. Reuses signature
# scores from 13_survival_signatures.py. See REFERENCES.md.
import os, numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; SV=os.path.join(ROOT,"Survival"); OUT=os.path.join(SV,"results")

sig=pd.read_csv(os.path.join(OUT,"tcga_signature_scores.csv"), index_col=0)[["TLS","JAK_STAT","TGFb"]]
surv=pd.read_csv(os.path.join(SV,"TCGA_BRCA","BRCA_survival.txt"), sep="\t").set_index("sample")
clin=pd.read_csv(os.path.join(SV,"TCGA_BRCA","BRCA_clinicalMatrix"), sep="\t", index_col=0)
tnbc=set(clin.index[(clin["ER_Status_nature2012"]=="Negative")&(clin["PR_Status_nature2012"]=="Negative")&(clin["HER2_Final_Status_nature2012"]=="Negative")])

df=sig.join(surv[["DFI","DFI.time","PFI","PFI.time"]])
endpoints={"PFI":("PFI","PFI.time"),"DFI":("DFI","DFI.time")}
groups={"ALL":set(df.index),"TNBC":tnbc}
res=[]
for ep,(ev,tm) in endpoints.items():
    for gname,gset in groups.items():
        for name in ["TLS","JAK_STAT","TGFb"]:
            d=df.loc[[i for i in df.index if i in gset],[name,ev,tm]].dropna()
            d=d[d[tm]>0].copy(); d["t"]=d[tm]/30.44; d["e"]=d[ev].astype(int)
            nev=int(d["e"].sum())
            if len(d)<20 or nev<5:
                res.append(dict(endpoint=ep,group=gname,signature=name,n=len(d),events=nev,note="underpowered")); continue
            d["z"]=(d[name]-d[name].mean())/d[name].std()
            cph=CoxPHFitter().fit(d[["z","t","e"]],"t","e")
            hr=float(np.exp(cph.params_["z"])); lo,hi=np.exp(cph.confidence_intervals_.loc["z"].values); p=float(cph.summary.loc["z","p"])
            grp=d[name]>=d[name].median()
            lr=logrank_test(d.loc[grp,"t"],d.loc[~grp,"t"],d.loc[grp,"e"],d.loc[~grp,"e"])
            res.append(dict(endpoint=ep,group=gname,signature=name,n=len(d),events=nev,HR_per_SD=round(hr,3),CI_low=round(lo,3),CI_high=round(hi,3),cox_p=round(p,4),logrank_p=round(lr.p_value,4)))
            if ep=="PFI":
                fig,ax=plt.subplots(figsize=(6,5)); km=KaplanMeierFitter()
                km.fit(d.loc[grp,"t"],d.loc[grp,"e"],label=f"{name} high (n={int(grp.sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
                km.fit(d.loc[~grp,"t"],d.loc[~grp,"e"],label=f"{name} low (n={int((~grp).sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
                ax.set_title(f"TCGA-BRCA {gname} — PFI by {name}\nHR/SD={hr:.2f} p={p:.3g}, log-rank {lr.p_value:.3g}"); ax.set_ylim(0,1); ax.set_xlabel("months")
                fig.tight_layout(); fig.savefig(os.path.join(OUT,f"KM_PFI_{gname}_{name}.png"),dpi=150); plt.close(fig)
out=pd.DataFrame(res); out.to_csv(os.path.join(OUT,"survival_recurrence_stats.csv"),index=False)
print(out.to_string(index=False)); print("RECURRENCE_DONE")
