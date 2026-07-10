# ============================================================================
# 03_label_transfer.R  —  Step 3: transfer cell-type labels onto the Visium spots
#
# Refactor of Label_transfer/label_transfer_script(step3).R.
# For each processed Visium sample, FindTransferAnchors + TransferData project
# the reference celltype_major labels into a "predictions" assay.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

ref <- readRDS(REF_RDS)

run_label_transfer <- function(spatial_obj, sample_name) {
  message("  label transfer: ", sample_name)
  anchors <- FindTransferAnchors(
    reference = ref, query = spatial_obj,
    normalization.method = "SCT", verbose = FALSE
  )
  predictions <- TransferData(
    anchorset = anchors, refdata = ref$celltype_major,
    prediction.assay = TRUE, weight.reduction = spatial_obj[["pca"]], dims = 1:30
  )
  spatial_obj[["predictions"]] <- predictions
  spatial_obj
}

for (id in SAMPLES) {
  message("=== Step 3 | ", id, " ===")
  obj <- readRDS(processed_rds(id))
  obj <- run_label_transfer(obj, id)
  saveRDS(obj, labeled_rds(id))
  message("  saved: ", labeled_rds(id))
}

message("Step 3 complete.")
