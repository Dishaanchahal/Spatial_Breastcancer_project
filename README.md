# Spatial Transcriptomics of Breast Cancer (Visium + Wu et al. 2021 reference)

Cell-type deconvolution of six 10x Visium breast-cancer slides by label transfer
from the Wu et al. 2021 scRNA-seq atlas (GSE176078), validated against
pathologist H&E annotations and summarized as tumor-microenvironment (TME)
composition per clinical subtype.

## Samples

| # | Sample ID | Subtype | PCA dims | Cluster res | Notes |
|---|-----------|---------|----------|-------------|-------|
| 1 | 1142243F  | TNBC | 30 | 0.8 | — |
| 2 | 1160920F  | TNBC | 30 | 0.8 | empty cluster 8 removed |
| 3 | CID4290   | ER+  | 20 | 0.6 | — |
| 4 | CID4465   | TNBC | 20 | 0.5 | — |
| 5 | CID44971  | TNBC | 20 | 0.5 | — |
| 6 | CID4535   | ER+  | 20 | 0.5 | — |

Reference cell types (`celltype_major`): Cancer Epithelial, T-cells, Myeloid,
CAFs, Endothelial, B-cells, Plasmablasts, PVL, Normal Epithelial.

## Pipeline

Run in order; each script sources `00_config.R` for all paths and parameters.

| Script | Step | Input | Output |
|--------|------|-------|--------|
| `00_config.R` | shared config (paths, sample list, per-sample params) | — | — |
| `01_process_spatial.R` | QC → SCTransform → cluster → spatially variable features | `samples/<ID>/` | `Processed_objects_step1/sample_<i>/visium_sample<i>.rds` |
| `02_process_reference.R` | build SCT-normalized scRNA-seq reference | `Wu_etal_2021_BRCA_scRNASeq/` | `ref_processed.rds` |
| `03_label_transfer.R` | FindTransferAnchors + TransferData | processed `.rds` + reference | `Label_transfer/sample<i>_<ID>_labeled.rds` |
| `04_label_transfer_plots.R` | spatial plots of 9 cell-type scores | labeled `.rds` | `Label_transfer/label_transfer_s<i>.png` |
| `05_pathologist_annotations.R` | overlay pathologist `Classification` | labeled `.rds` + `metadata/` | `Pathologists_annotations/…png` |
| `06_validate_annotations.R` | mean prediction score per pathologist class (heatmaps) | labeled `.rds` + `metadata/` | `Validating/validation_heatmap_s<i>.png` |
| `07_tme_comparison.R` | TME composition, computed from objects | labeled `.rds` | `TME_comparison/tme_comparison*.png` + `tme_summary_computed.csv` |
| `03b_deconvolution_rctd.R` | **proper** spot deconvolution (RCTD) — see notes | processed `.rds` + reference | `Deconvolution_RCTD/weights_<ID>.csv` |

## Running on Explorer

```bash
cd /scratch/chahal.d/Spatial_GenAI
sbatch run_pipeline.sh          # full pipeline as a SLURM job
# or a single step:
module load R/4.4.1
Rscript 07_tme_comparison.R
```

R packages required: `Seurat`, `ggplot2`, `patchwork`, `dplyr`, `reshape2`,
`pheatmap`, `scales`, `RColorBrewer`.

## What changed vs. the original scripts

The original exploratory scripts are preserved unchanged in
`_original_scripts/`. This refactor fixes the reproducibility issues:

1. **Consistent paths.** All I/O goes through `00_config.R`. The original scripts
   saved step-1 output to the working directory but read it from
   `Processed_objects_spatial/…` — paths that did not match the real folders.
2. **Renamed folder.** `Processed_objects_spatial--- step1` (spaces + dashes)
   → `Processed_objects_step1`.
3. **De-duplicated step 1.** Six copy-pasted sample blocks → one loop over
   `SAMPLES`/`PARAMS`; each sample's original tuning is preserved.
4. **Robust cluster removal.** Idents-based subsetting replaces the fragile
   `seurat_clusters != 8` factor comparison (the "cluster 8 removal failed" note).
5. **Self-contained validation.** `06_validate_annotations.R` loads its own
   objects/metadata instead of relying on variables left in the session.
6. **TME computed, not hardcoded.** `07_tme_comparison.R` derives the per-sample
   cell-type proportions from the labeled objects (was a typed-in table).
7. **Spot QC added.** Step 1 now computes mitochondrial % and filters spots
   (min genes/counts, max mito %) uniformly, and regresses mito in SCTransform.
   The original did no spot QC.

## Methodological notes (literature-validated)

- **Deconvolution vs. label transfer.** Visium spots (55 um) are multicellular,
  so each spot is a *mixture* of cell types. Seurat anchor label transfer
  (`TransferData`) produces per-spot prediction *scores* — a useful quick view,
  but independent benchmarks show dedicated deconvolution methods
  (cell2location, RCTD, CARD, SpatialDWLS/Tangram) are the appropriate tools for
  estimating spot composition, consistently topping accuracy benchmarks across
  many datasets (Li et al. 2022, *Nat Methods*; Li et al. 2023, *Nat Commun*;
  Sang-aram et al. 2024, *eLife*). `03b_deconvolution_rctd.R` adds RCTD as the
  corrected composition estimate; treat the step-7 TME figure as an
  approximation until RCTD weights are available.
- **Biological plausibility.** The subtype contrast in the data — ER+ samples
  dominated by cancer epithelium with sparse immune content, TNBC samples more
  immune/stroma-infiltrated (T-cells, CAFs) — is consistent with reports that
  TNBC carries higher immune infiltration (CD8+ T cells alongside
  immunosuppressive Tregs/M2 macrophages) and stromal remodeling than
  hormone-receptor-positive breast cancer.

## Suggested next steps

1. Install `spacexr` and run `03b_deconvolution_rctd.R`; rebuild the TME
   comparison from RCTD weights for a defensible composition estimate.
2. Add statistics to the TME contrast (per-cell-type proportions, TNBC vs ER+)
   rather than eyeballing stacked bars — with n=4 vs 2 this is descriptive only,
   so state that explicitly.
3. Spatial neighborhood / colocalization analysis (e.g. are T-cells adjacent to
   tumor epithelium vs excluded?) to address immune exclusion, the key
   subtype-distinguishing feature in the TNBC TME literature.
4. Harmonize QC/clustering parameters across samples (or justify per-sample
   choices) and record package versions (`sessionInfo()` / renv lockfile).

## RCTD cross-validation (Option 3)

An independent deconvolution with RCTD (`spacexr`) was run as a method cross-check
(`08_rctd_validation.R` → `09_compare_c2l_rctd.R` → `cell2location/09b_plot_concordance.py`).

- Per-sample composition from **cell2location and RCTD agree strongly: Pearson r = 0.93**.
- Systematic differences: RCTD assigns **more Cancer Epithelial and CAFs**, cell2location
  assigns **more Normal Epithelial** — consistent with the earlier caveat that cell2location
  may over-call "Normal Epithelial" within tumor regions (RCTD reallocates it to Cancer Epithelial).
- Subtype-level conclusions hold under **both** methods: ER+ tumor-dominated/immune-cold;
  TNBC more CAF/immune-infiltrated (CID4465 the most stroma/immune-rich).

See `Validation_RCTD/c2l_vs_rctd_scatter.png` and `Validation_RCTD/rctd_composition.csv`.

## Cell–cell communication (Option 2)

Ligand–receptor communication among the 9 major cell types was inferred with
**LIANA** (consensus rank-aggregate) on the scRNA-seq reference
(`CellComm/10_cell_communication.py`). 15,616 interactions ranked.

Key findings (see `CellComm/liana_dotplot_structural_to_immune.png`,
`CellComm/liana_key_chemokines.csv`):
- **CXCL12 → CXCR4** (stroma/endothelium/myeloid → T- and B-cells) is prominent and,
  for **CAFs → T/B cells, highly specific** — the canonical immune recruitment/retention
  axis, matching the CAF-immune interplay seen in the CAF-high TNBC sample (CID4465).
- **Antigen presentation**: MHC-II → CD4 (HLA-DP/DR → CD4) and MHC-I → CD8A from CAFs,
  PVL, and epithelium to T-cells.
- **ECM signaling**: COL1A1/COL1A2 → CD44 and VIM → CD44 from CAFs/stroma to immune cells.
- **B-lineage support**: TNFSF13B (BAFF) → Plasmablasts/T-cells from CAFs and epithelium.

## Tissue architecture — spatial niches + TLS (Option: architecture)

Spatial niches were derived by clustering neighborhood-aggregated cell2location
composition (k-NN, k=15; KMeans, 8 niches) — the cellular-neighborhood approach
(Goltsev 2018; Schürch 2020; Janesick 2023) implemented with Squidpy-style logic (Palla 2022).
Script: `TissueArchitecture/11_tissue_architecture.py`.

- A **lymphoid niche** (niche 4: B-cell + T-cell rich) was identified whose location
  coincides with peaks of the **Cabrita 2020** 9-gene **TLS signature** — a putative
  tertiary lymphoid structure, most prominent in **CID4465 (TNBC)**. Two independent
  methods (composition niche + TLS gene signature) agree on its location.
- **Immune exclusion** (`immune_exclusion_and_cxcl12_cxcr4.csv`): T-cells co-localize with
  **CAFs in every sample** (Spearman 0.25–0.71); the spatial CXCL12(neighborhood)→CXCR4
  coupling is modest-positive, consistent with the LIANA CAF→T/B recruitment axis.

## Pathway activity — PROGENy (Option: pathway)

PROGENy pathway activity (MLM) via decoupleR — `TissueArchitecture/12_pathway_activity.py`
(Schubert 2018; Badia-i-Mompel 2022).

- **By niche:** JAK-STAT is highest in the **lymphoid/TLS niche** (interferon/cytokine);
  TGFβ is highest in the **CAF-rich tumor-stroma niche** (fibrosis) — pathway activity maps
  cleanly onto tissue architecture.
- **By subtype:** TNBC shows higher **JAK-STAT** than ER+, consistent with a more inflamed TME.
  Descriptive only (n = 4 TNBC vs 2 ER+).

All method and biology citations are in **REFERENCES.md**.

## External cross-validation (independent 10x dataset)

The pipeline was applied to an **independent** 10x Visium breast-cancer dataset
(Block A, Sections 1–2; invasive ductal carcinoma; 3,798 + 3,987 spots; 18,036 genes
shared with our reference), **reusing our reference signatures without retraining**.
Scripts: `ExternalValidation/` (`xval_c2l.py`, `xval_downstream.py`).

- **Composition reproduces:** same 9-cell-type repertoire, tumor-dominated
  (Cancer Epithelial ~0.55), and the two adjacent sections are near-identical
  (internal consistency).
- **Pathway structure reproduces:** JAK-STAT strongly active (~12), and the full
  PROGENy directional signature matches our cohort (positive JAK-STAT/Hypoxia/TGFβ/TNFα;
  negative p53/PI3K/MAPK/Estrogen).
- **TLS detection generalizes:** focal Cabrita-signature TLS foci are detected
  (consistent between sections), not diffuse.

Conclusion: the deconvolution reference, TLS signature, and pathway analysis
transfer to data the pipeline was not built on. See REFERENCES.md.

## Clinical validation — survival (TCGA-BRCA)

The TLS (Cabrita 2020) and PROGENy **JAK-STAT** / **TGFβ** signatures — the axes we
found spatially (TLS/lymphoid niche; JAK-STAT in the lymphoid niche; TGFβ in the CAF niche) —
were projected onto **TCGA-BRCA** bulk RNA-seq (UCSC Xena; n = 1,201) and tested for
overall survival (Cox HR per SD; KM median split). Script: `Survival/13_survival_signatures.py`.

| Signature | HR per SD (95% CI) | Cox p | direction |
|---|---|---|---|
| TLS       | 0.89 (0.78–1.02) | 0.10  | higher → better OS (trend) |
| JAK-STAT  | 0.88 (0.76–1.01) | 0.065 | higher → better OS (borderline) |
| TGFβ      | 1.00 (0.87–1.16) | 0.96  | neutral |

**Honest interpretation:** the immune/TLS signatures trend toward *better* survival
(HR < 1), in the expected direction, but **none reach statistical significance** in
TCGA-BRCA overall survival — expected, since TCGA-BRCA OS is underpowered (few death
events, mixed subtypes). The JAK-STAT-high curve separates favorably for ~10 years
before a small-number tail. Stronger tests: **METABRIC** (larger, longer follow-up),
subtype-stratified (basal/TNBC), or recurrence endpoints (DFI/PFI). See REFERENCES.md.

### Subtype-stratified survival (TCGA-BRCA)

Rerun within **Basal** (PAM50; n=140, 20 events) and **TNBC** (ER-/PR-/HER2-; n=123, 19 events).
Script: `Survival/14_survival_subtype.py`.

| Subset | Signature | HR per SD | Cox p |
|---|---|---|---|
| TNBC  | JAK-STAT | 0.73 | 0.19 |
| TNBC  | TLS      | 0.81 | 0.38 |
| Basal | JAK-STAT | 0.96 | 0.85 |
| Basal | TLS      | 0.90 | 0.63 |

In **TNBC the protective effect sizes are larger** than the full cohort (JAK-STAT HR 0.73,
TLS 0.81 vs 0.88/0.89) — directionally consistent with TME/TLS mattering more in TNBC — but
with only ~19 death events the TNBC subset is **severely underpowered and non-significant**.
Basal near-null. Honest conclusion: effect direction moves as hypothesized, but TCGA subtype
subsets lack events to confirm; needs METABRIC / pooled cohorts / recurrence (DFI/PFI) endpoints.

### Recurrence-endpoint survival (TCGA-BRCA) — significant

Overall survival is event-poor in TCGA-BRCA; using **recurrence endpoints** (PFI/DFI, more
events) the signatures reach significance (continuous Cox HR per SD). Script:
`Survival/15_survival_recurrence.py`.

| Endpoint | Signature | HR per SD (95% CI) | Cox p |
|---|---|---|---|
| PFI | **TLS**  | 0.82 (0.71–0.96) | **0.011** |
| DFI | **TLS**  | 0.78 (0.63–0.95) | **0.016** |
| DFI | **TGFβ** | 1.26 (1.02–1.56) | **0.034** |
| PFI | TGFβ | 1.15 (0.98–1.35) | 0.080 |
| PFI | JAK-STAT | 0.87 (0.74–1.01) | 0.068 |

**Clinically validates the spatial findings:** higher **TLS** (the lymphoid niche) → significantly
**less progression/recurrence** (protective); higher **TGFβ** (the CAF niche) → significantly
**more recurrence** (harmful). JAK-STAT protective trend. Continuous Cox is the primary test
(median-split log-rank is less powerful). TNBC subset is directionally consistent (JAK-STAT
HR~0.72–0.78) but underpowered (~16–23 events).

### METABRIC validation (independent cohort) — TLS replicates strongly

Independent, well-powered cohort (n=1,903; **1,103 death events**; OS). Script:
`Survival/16_survival_metabric.py`.

| Group | Signature | HR per SD (95% CI) | Cox p |
|---|---|---|---|
| ALL | **TLS** | 0.87 (0.82–0.93) | **3e-5** |
| ALL | JAK-STAT | 1.04 (0.98–1.10) | 0.20 (ns) |
| ALL | TGFβ | 1.05 (0.99–1.11) | 0.14 (ns) |

**Cross-cohort verdict:**
- **TLS is robustly, significantly protective** — METABRIC OS (p=3e-5) *replicates* the TCGA
  recurrence result (PFI p=0.011, DFI p=0.016). Two independent cohorts, two endpoint types →
  the TLS/lymphoid-niche prognostic effect is real and reproducible.
- **JAK-STAT did NOT replicate** (METABRIC HR 1.04, ns; TCGA had only a borderline trend) —
  not a robust prognostic signal; dropped as a claim.
- **TGFβ** is harmful in direction in both cohorts but significant only in TCGA DFI (HR 1.26) —
  suggestive, not confirmed.
Subtype subsets (Basal, ER-/HER2-) show TLS in the same protective direction (HR~0.88) but are
underpowered within-subset.
