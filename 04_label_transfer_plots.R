# ============================================================================
# 04_label_transfer_plots.R  —  Step 4: spatial plots of cell-type predictions
#
# Refactor of Label_transfer/label_transfer_plot_script(step4).R.
# Plots all 9 cell-type prediction scores spatially for each sample.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

for (id in SAMPLES) {
  i <- match(id, SAMPLES)
  message("=== Step 4 | ", id, " ===")

  obj <- readRDS(labeled_rds(id))
  DefaultAssay(obj) <- "predictions"

  p <- SpatialFeaturePlot(
        obj, features = CELL_TYPES, ncol = 3,
        alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 0.5
      ) & scale_fill_gradientn(
        colors = c("navy", "blue", "cyan", "yellow", "red"), limits = c(0, 0.5)
      )

  ggsave(file.path(DIR_LABELED, paste0("label_transfer_s", i, ".png")),
         plot = p, width = 14, height = 14, dpi = 300)
}

message("Step 4 complete.")
