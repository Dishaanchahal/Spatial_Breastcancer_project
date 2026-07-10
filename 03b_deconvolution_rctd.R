# ============================================================================
# 03b_deconvolution_rctd.R  —  Proper spot deconvolution with RCTD (spacexr)
#
# WHY THIS EXISTS (literature-validated correction):
# 10x Visium spots are 55 um and contain multiple cells, so a spot's expression
# is a MIXTURE of cell types. Seurat anchor label transfer (steps 3/4/7) gives
# per-spot prediction *scores*, which approximate but are not a validated
# deconvolution. Independent benchmarks recommend dedicated deconvolution
# methods (cell2location, RCTD, CARD, SpatialDWLS) for estimating spot
# composition (Li et al. 2022 Nat Methods; Li et al. 2023 Nat Commun;
# Sang-aram et al. 2024 eLife). RCTD is a top performer and is pure R.
#
# REQUIRES: the 'spacexr' package (github.com/dmcable/spacexr). Not yet
# installed on Explorer — install into your R library first:
#   module load R/4.4.1
#   R -e 'remotes::install_github("dmcable/spacexr")'
#
# NOTE: coordinate/slot accessors differ slightly across Seurat versions;
# validate the GetTissueCoordinates() output columns on first run.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(spacexr) })

# ---- Build the RCTD reference from the scRNA-seq atlas (once) ----
ref <- readRDS(REF_RDS)
ref_counts <- GetAssayData(ref, assay = "RNA", slot = "counts")
cell_types <- as.factor(ref$celltype_major); names(cell_types) <- colnames(ref)
nUMI_ref   <- ref$nCount_RNA;                names(nUMI_ref)   <- colnames(ref)
reference  <- Reference(ref_counts, cell_types, nUMI_ref)

for (id in SAMPLES) {
  i <- match(id, SAMPLES)
  message("=== Step 3b (RCTD) | ", id, " ===")

  obj    <- readRDS(processed_rds(id))
  counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")

  # tissue coordinates -> data.frame(x, y) with barcode rownames
  coords <- as.data.frame(GetTissueCoordinates(obj))
  coords <- coords[, sapply(coords, is.numeric)][, 1:2]
  colnames(coords) <- c("x", "y")

  puck <- SpatialRNA(coords, counts, colSums(counts))

  myRCTD <- create.RCTD(puck, reference, max_cores = 4)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

  # normalized per-spot cell-type weights (spots x cell types)
  weights <- as.matrix(myRCTD@results$weights)
  weights <- sweep(weights, 1, rowSums(weights), "/")

  saveRDS(myRCTD, file.path(DIR_DECONV, paste0("rctd_", id, ".rds")))
  write.csv(weights, file.path(DIR_DECONV, paste0("weights_", id, ".csv")))
  message("  saved RCTD weights: ", nrow(weights), " spots x ", ncol(weights), " types")
}

message("Step 3b complete. Use these weights for TME composition (preferred over prediction scores).")
