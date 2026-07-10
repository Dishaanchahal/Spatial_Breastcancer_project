# ============================================================================
# 00_config.R  —  Central configuration for the spatial transcriptomics pipeline
# Breast cancer Visium + Wu et al. 2021 scRNA-seq reference (GSE176078)
#
# Every pipeline script sources this file, so all paths, sample IDs and
# per-sample parameters live in ONE place. Edit here, not in the scripts.
# ============================================================================

PROJECT_ROOT <- "/scratch/chahal.d/Spatial_GenAI"

# ---- Samples (fixed order) and their clinical subtype -----------------------
SAMPLES  <- c("1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535")
SUBTYPES <- c("TNBC",     "TNBC",     "ER+",     "TNBC",    "TNBC",     "ER+")
names(SUBTYPES) <- SAMPLES

# The 9 major cell types transferred from the reference (celltype_major)
CELL_TYPES <- c("Cancer Epithelial", "T-cells", "Myeloid", "CAFs", "Endothelial",
                "B-cells", "Plasmablasts", "PVL", "Normal Epithelial")

# ---- Spot-level QC thresholds (applied uniformly in step 1) -----------------
# The original scripts did no QC filtering (spots were only inspected visually
# and one cluster was dropped by hand). These are sensible starting defaults for
# tumor Visium; tune per tissue. Mito genes are symbols ("MT-*") in this data.
QC <- list(min_features = 200, min_counts = 500, max_mt_percent = 25)

# ---- Per-sample step-1 parameters (preserved from the original hand-tuning) --
# dims        : PCA dims used for neighbors / clusters / UMAP
# resolution  : FindClusters resolution
# remove_clusters : QC-driven clusters to drop before downstream analysis
PARAMS <- list(
  "1142243F" = list(dims = 30, resolution = 0.8, remove_clusters = integer(0)),
  "1160920F" = list(dims = 30, resolution = 0.8, remove_clusters = c(8)),  # empty cluster 8
  "CID4290"  = list(dims = 20, resolution = 0.6, remove_clusters = integer(0)),
  "CID4465"  = list(dims = 20, resolution = 0.5, remove_clusters = integer(0)),
  "CID44971" = list(dims = 20, resolution = 0.5, remove_clusters = integer(0)),
  "CID4535"  = list(dims = 20, resolution = 0.5, remove_clusters = integer(0))
)

# ---- Directories ------------------------------------------------------------
DIR_SAMPLES  <- file.path(PROJECT_ROOT, "samples")                       # raw Visium inputs
DIR_METADATA <- file.path(PROJECT_ROOT, "metadata")                      # pathologist annotations
DIR_REF      <- file.path(PROJECT_ROOT, "Wu_etal_2021_BRCA_scRNASeq")    # scRNA-seq reference
DIR_STEP1    <- file.path(PROJECT_ROOT, "Processed_objects_step1")       # per-sample processed objects
DIR_LABELED  <- file.path(PROJECT_ROOT, "Label_transfer")               # label-transferred objects + plots
DIR_ANNOT    <- file.path(PROJECT_ROOT, "Pathologists_annotations")
DIR_VALID    <- file.path(PROJECT_ROOT, "Validating")
DIR_TME      <- file.path(PROJECT_ROOT, "TME_comparison")
DIR_DECONV   <- file.path(PROJECT_ROOT, "Deconvolution_RCTD")   # proper deconvolution (step 3b)

REF_RDS      <- file.path(DIR_REF, "ref_processed.rds")

# ---- Path helpers -----------------------------------------------------------
processed_rds <- function(id) {
  i <- match(id, SAMPLES)
  file.path(DIR_STEP1, paste0("sample_", i), paste0("visium_sample", i, ".rds"))
}

labeled_rds <- function(id) {
  i <- match(id, SAMPLES)
  file.path(DIR_LABELED, paste0("sample", i, "_", id, "_labeled.rds"))
}

sample_metadata <- function(id) file.path(DIR_METADATA, paste0(id, "_metadata.csv"))

# Create output dirs if missing (safe to call repeatedly)
for (d in c(DIR_STEP1, DIR_LABELED, DIR_ANNOT, DIR_VALID, DIR_TME, DIR_DECONV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
