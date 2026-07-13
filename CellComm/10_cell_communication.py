# Cell-cell communication among the 9 major cell types (LIANA consensus,
# run on the scRNA-seq reference where per-cell expression is clean).
import os, scanpy as sc, liana as li
ROOT = "/scratch/chahal.d/Spatial_GenAI"; OUT = os.path.join(ROOT, "CellComm"); os.makedirs(OUT, exist_ok=True)

adata = sc.read_h5ad(os.path.join(ROOT, "cell2location", "data", "ref.h5ad"))
adata.X = adata.X.astype("float32")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("reference:", adata.shape, "| cell types:", sorted(adata.obs["celltype_major"].unique()))

li.mt.rank_aggregate(adata, groupby="celltype_major", expr_prop=0.1, use_raw=False, verbose=True)
res = adata.uns["liana_res"]
res.to_csv(os.path.join(OUT, "liana_rank_aggregate.csv"), index=False)
print("total ranked interactions:", len(res))

# key recruitment/immune chemokine axes flagged by the literature
import re
pat = r"CXCL12|CXCR4|CXCL9|CXCL10|CXCL11|CXCR3|IL15|CCL19|CCL5|CXCL13|CCL21"
mask = res["ligand_complex"].str.contains(pat, case=False, na=False) | \
       res["receptor_complex"].str.contains(pat, case=False, na=False)
res[mask].to_csv(os.path.join(OUT, "liana_key_chemokines.csv"), index=False)

# dotplot: structural/myeloid sources -> lymphoid targets
struct = ["CAFs", "Endothelial", "PVL", "Myeloid", "Cancer Epithelial", "Normal Epithelial"]
immune = ["T-cells", "B-cells", "Plasmablasts", "Myeloid"]
try:
    p = li.pl.dotplot(adata, colour="magnitude_rank", size="specificity_rank",
                      source_labels=struct, target_labels=immune, top_n=25,
                      orderby="magnitude_rank", orderby_ascending=True, figure_size=(12, 11))
    p.save(os.path.join(OUT, "liana_dotplot_structural_to_immune.png"), dpi=150, verbose=False)
    print("wrote dotplot")
except Exception as e:
    print("dotplot step failed (results CSV still written):", repr(e))
print("CCC_DONE")
