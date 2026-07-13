# Tissue architecture: spatial niches, TLS scoring, immune exclusion, and
# spatial CXCL12-CXCR4 colocalization.
# Methods/refs (see REFERENCES.md): cellular-neighborhood clustering
# (Goltsev 2018; Janesick 2023; Oliveira 2025), Squidpy (Palla 2022),
# TLS 9-gene signature (Cabrita 2020; Helmink 2020), COMMOT concept (Cang 2023).
import os, numpy as np, pandas as pd, scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/scratch/chahal.d/Spatial_GenAI"; OUT=os.path.join(ROOT,"TissueArchitecture"); os.makedirs(OUT,exist_ok=True)
CELL_TYPES=["Cancer Epithelial","T-cells","Myeloid","CAFs","Endothelial","B-cells","Plasmablasts","PVL","Normal Epithelial"]
TLS_GENES=["CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS"]  # Cabrita 2020 Nature
K=15; N_NICHES=8

# --- load cell2location abundances (obs) + coords + sample ---
vm=sc.read_h5ad(os.path.join(ROOT,"cell2location","results","vis_mapped.h5ad"))
ab=vm.obs[CELL_TYPES].to_numpy(float)
coords=np.asarray(vm.obsm["spatial"],float); samples=vm.obs["sample"].astype(str).to_numpy()
sub=dict(zip(vm.obs["sample"].astype(str), vm.obs["subtype"].astype(str)))

# --- neighborhood-aggregated composition, then KMeans -> niches ---
nbr_comp=np.zeros_like(ab)
for s in pd.unique(samples):
    m=np.where(samples==s)[0]
    nn=NearestNeighbors(n_neighbors=min(K,len(m))).fit(coords[m])
    _,idx=nn.kneighbors(coords[m])
    nbr_comp[m]=ab[m][idx].mean(1)
km=KMeans(n_clusters=N_NICHES,random_state=0,n_init=10).fit(nbr_comp)
niche=km.labels_.astype(int)
pd.DataFrame({"barcode":vm.obs_names,"sample":samples,"niche":niche}).to_csv(os.path.join(OUT,"niche_labels.csv"),index=False)

# niche composition heatmap (mean own-spot abundance per niche)
comp=pd.DataFrame(ab,columns=CELL_TYPES); comp["niche"]=niche
nmean=comp.groupby("niche").mean()
nmean.to_csv(os.path.join(OUT,"niche_composition.csv"))
fig,ax=plt.subplots(figsize=(9,6)); im=ax.imshow(nmean.values,aspect="auto",cmap="magma")
ax.set_xticks(range(len(CELL_TYPES))); ax.set_xticklabels(CELL_TYPES,rotation=45,ha="right")
ax.set_yticks(range(N_NICHES)); ax.set_yticklabels([f"niche {i}" for i in range(N_NICHES)])
ax.set_title("Spatial niche composition (mean cell2location abundance)"); fig.colorbar(im)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"niche_composition_heatmap.png"),dpi=150); plt.close(fig)

# niche spatial maps
for s in pd.unique(samples):
    m=np.where(samples==s)[0]; xy=coords[m]
    fig,ax=plt.subplots(figsize=(6,6))
    sctr=ax.scatter(xy[:,0],-xy[:,1],c=niche[m],cmap="tab10",s=8,vmin=0,vmax=9)
    ax.set_title(f"{s} ({sub[s]}) niches"); ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
    fig.colorbar(sctr,label="niche"); fig.tight_layout()
    fig.savefig(os.path.join(OUT,f"niches_{s}.png"),dpi=140); plt.close(fig)

# --- TLS score (Cabrita signature) from full-gene expression ---
ge=sc.read_h5ad(os.path.join(ROOT,"cell2location","data","vis.h5ad"))
ge=ge[vm.obs_names].copy()
sc.pp.normalize_total(ge,target_sum=1e4); sc.pp.log1p(ge)
present=[g for g in TLS_GENES if g in ge.var_names]
sc.tl.score_genes(ge,present,score_name="TLS")
tls=ge.obs["TLS"].to_numpy()
print("TLS genes present:",present)
for s in pd.unique(samples):
    m=np.where(samples==s)[0]; xy=coords[m]
    fig,ax=plt.subplots(figsize=(6,6))
    sctr=ax.scatter(xy[:,0],-xy[:,1],c=tls[m],cmap="viridis",s=8)
    ax.set_title(f"{s} ({sub[s]}) TLS score"); ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect("equal")
    fig.colorbar(sctr); fig.tight_layout(); fig.savefig(os.path.join(OUT,f"tls_{s}.png"),dpi=140); plt.close(fig)
tdf=pd.DataFrame({"niche":niche,"TLS":tls}).groupby("niche").mean()
tdf.to_csv(os.path.join(OUT,"tls_by_niche.csv"))

# --- immune exclusion + CXCL12-CXCR4 spatial colocalization (per sample) ---
def gene_vec(g): return np.asarray(ge[:,g].X.todense()).ravel() if g in ge.var_names else None
cxcl12=gene_vec("CXCL12"); cxcr4=gene_vec("CXCR4")
rows=[]
for s in pd.unique(samples):
    m=np.where(samples==s)[0]
    T=ab[m,CELL_TYPES.index("T-cells")]; CE=ab[m,CELL_TYPES.index("Cancer Epithelial")]; CF=ab[m,CELL_TYPES.index("CAFs")]
    r_ce=spearmanr(T,CE).correlation; r_cf=spearmanr(T,CF).correlation
    # distance of each spot to nearest tumor spot -> corr with T-cell abundance
    tum=np.where(CE>=np.quantile(CE,0.75))[0]
    d=NearestNeighbors(n_neighbors=1).fit(coords[m][tum]); dist,_=d.kneighbors(coords[m]); dist=dist.ravel()
    r_dist=spearmanr(T,dist).correlation
    # CXCL12(neighbor-mean) vs CXCR4(own)
    r_lr=np.nan
    if cxcl12 is not None and cxcr4 is not None:
        nn=NearestNeighbors(n_neighbors=min(K,len(m))).fit(coords[m]); _,idx=nn.kneighbors(coords[m])
        lig_nb=cxcl12[m][idx].mean(1)
        r_lr=pearsonr(lig_nb,cxcr4[m])[0]
    rows.append(dict(sample=s,subtype=sub[s],
        spearman_T_vs_CancerEpi=r_ce, spearman_T_vs_CAF=r_cf,
        spearman_T_vs_dist_to_tumor=r_dist, pearson_CXCL12nb_vs_CXCR4=r_lr))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"immune_exclusion_and_cxcl12_cxcr4.csv"),index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False))
print("TISSUE_ARCH_DONE")
