# Subtype-stratified survival on TCGA-BRCA: rerun within Basal (PAM50) and TNBC
# (ER-/PR-/HER2-), where TME/TLS biology is most relevant. Reuses the signature
# scores computed in 13_survival_signatures.py. See REFERENCES.md.
import os, numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/scratch/chahal.d/Spatial_GenAI"; SV=os.path.join(ROOT,"Survival"); OUT=os.path.join(SV,"results")

sig=pd.read_csv(os.path.join(OUT,"tcga_signature_scores.csv"), index_col=0)   # TLS,JAK_STAT,TGFb,OS,OS.time,OS_months
clin=pd.read_csv(os.path.join(SV,"TCGA_BRCA","BRCA_clinicalMatrix"), sep="\t", index_col=0)
pam=clin["PAM50Call_RNAseq"]; er=clin["ER_Status_nature2012"]; pr=clin["PR_Status_nature2012"]; her2=clin["HER2_Final_Status_nature2012"]

df=sig.copy()
df["Basal"]=[pam.get(s,"")=="Basal" for s in df.index]
df["TNBC"] =[(er.get(s,"")=="Negative") and (pr.get(s,"")=="Negative") and (her2.get(s,"")=="Negative") for s in df.index]
print("n total=%d | Basal=%d | TNBC=%d" % (len(df), df["Basal"].sum(), df["TNBC"].sum()))

subsets={"ALL":df.index, "Basal":df.index[df["Basal"]], "TNBC":df.index[df["TNBC"]]}
res=[]
for sub,idx in subsets.items():
    base=df.loc[idx]
    for name in ["TLS","JAK_STAT","TGFb"]:
        d=base[[name,"OS","OS_months"]].dropna(); d=d[d["OS_months"]>0].copy()
        ev=int(d["OS"].sum())
        if len(d)<20 or ev<5:
            res.append(dict(subset=sub,signature=name,n=len(d),events=ev,note="underpowered")); continue
        d["z"]=(d[name]-d[name].mean())/d[name].std()
        cph=CoxPHFitter().fit(d[["z","OS_months","OS"]],"OS_months","OS")
        hr=float(np.exp(cph.params_["z"])); lo,hi=np.exp(cph.confidence_intervals_.loc["z"].values); p=float(cph.summary.loc["z","p"])
        grp=d[name]>=d[name].median()
        lr=logrank_test(d.loc[grp,"OS_months"],d.loc[~grp,"OS_months"],d.loc[grp,"OS"],d.loc[~grp,"OS"])
        res.append(dict(subset=sub,signature=name,n=len(d),events=ev,HR_per_SD=round(hr,3),CI_low=round(lo,3),CI_high=round(hi,3),cox_p=round(p,4),logrank_p=round(lr.p_value,4)))
        if sub!="ALL":
            fig,ax=plt.subplots(figsize=(6,5)); km=KaplanMeierFitter()
            km.fit(d.loc[grp,"OS_months"],d.loc[grp,"OS"],label=f"{name} high (n={int(grp.sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
            km.fit(d.loc[~grp,"OS_months"],d.loc[~grp,"OS"],label=f"{name} low (n={int((~grp).sum())})"); km.plot_survival_function(ax=ax,ci_show=False)
            ax.set_title(f"TCGA-BRCA {sub} — OS by {name}\nHR/SD={hr:.2f} (p={p:.3g}), log-rank p={lr.p_value:.3g}"); ax.set_ylim(0,1); ax.set_xlabel("months")
            fig.tight_layout(); fig.savefig(os.path.join(OUT,f"KM_{sub}_{name}.png"),dpi=150); plt.close(fig)
out=pd.DataFrame(res); out.to_csv(os.path.join(OUT,"survival_subtype_stats.csv"),index=False)
print(out.to_string(index=False)); print("SUBTYPE_SURVIVAL_DONE")
