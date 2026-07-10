# ============================================================================
# 05_pathologist_annotations.R  —  Overlay pathologist H&E annotations
#
# Refactor of Pathologists_annotations/Pathologists_annotations.R.
# A pathologist classified each spot from the H&E image (Classification column
# in metadata). Here we attach those labels to the labeled objects and plot
# them spatially, to visually compare against the transferred cell types.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

plots <- list()
for (id in SAMPLES) {
  message("=== Step 5 | ", id, " ===")
  meta <- read.csv(sample_metadata(id), row.names = 1)
  message("  pathologist Classification counts:")
  print(table(meta$Classification))

  obj <- readRDS(labeled_rds(id))
  obj <- AddMetaData(obj, meta)

  plots[[id]] <- SpatialDimPlot(obj, group.by = "Classification") +
                 ggtitle(paste0(id, " (", SUBTYPES[id], ")"))
}

combined <- wrap_plots(plots, ncol = 3)
ggsave(file.path(DIR_ANNOT, "pathologist_annotations_all.png"),
       plot = combined, width = 18, height = 12, dpi = 300)

message("Step 5 complete.")
