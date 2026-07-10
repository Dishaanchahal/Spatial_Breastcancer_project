# ============================================================================
# 06_validate_annotations.R  —  Quantitative validation of predictions
#
# Refactor of Validating/validating_annotations.R. Now SELF-CONTAINED: it loads
# the labeled objects and metadata itself (the original relied on variables left
# in the session by the pathologist-annotation script).
#
# For each sample it builds a heatmap of the mean cell-type prediction score
# within each pathologist Classification. Concordance (e.g. high T-/B-cell
# scores in "Lymphocytes" spots) validates the label transfer.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(pheatmap) })

for (id in SAMPLES) {
  i <- match(id, SAMPLES)
  message("=== Step 6 | ", id, " ===")

  obj  <- readRDS(labeled_rds(id))
  meta <- read.csv(sample_metadata(id), row.names = 1)

  # spots x cell-type prediction scores
  df <- as.data.frame(t(GetAssayData(obj, assay = "predictions")))
  df$max <- NULL  # TransferData adds a 'max' row; not a cell type
  df$Classification <- meta[rownames(df), "Classification"]
  df <- df[!is.na(df$Classification), , drop = FALSE]

  # mean prediction score per pathologist class
  agg <- aggregate(. ~ Classification, data = df, FUN = mean)
  rownames(agg) <- agg$Classification
  agg$Classification <- NULL

  pheatmap(
    agg, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
    color = colorRampPalette(c("navy", "white", "red"))(50),
    main = paste0("Mean prediction score - ", id, " (", SUBTYPES[id], ")"),
    filename = file.path(DIR_VALID, paste0("validation_heatmap_s", i, ".png")),
    width = 10, height = 6
  )
}

message("Step 6 complete.")
