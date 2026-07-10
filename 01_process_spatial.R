# ============================================================================
# 01_process_spatial.R  —  Step 1: per-sample Visium processing
#
# Refactor of the original first_sample_script.R: the six copy-pasted sample
# blocks are now a single loop driven by SAMPLES / PARAMS in 00_config.R.
# Each sample's original tuning (dims, resolution, cluster removal) is preserved.
# Outputs go to Processed_objects_step1/sample_<i>/.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})

for (id in SAMPLES) {
  i <- match(id, SAMPLES)
  p <- PARAMS[[id]]
  message("=== Step 1 | sample ", i, ": ", id, " (", SUBTYPES[id], ") ===")

  out_dir <- file.path(DIR_STEP1, paste0("sample_", i))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Load count matrix + tissue image ----
  counts <- Read10X(file.path(DIR_SAMPLES, id, "filtered_feature_bc_matrix"), gene.column = 1)
  visium <- CreateSeuratObject(counts, assay = "Spatial")
  image  <- Read10X_Image(file.path(DIR_SAMPLES, id, "spatial"))
  image  <- image[Cells(visium)]
  DefaultAssay(image) <- "Spatial"
  visium[[id]] <- image

  # ---- Mitochondrial content + spot-level QC filtering ----
  # (Added: the original pipeline did no spot QC.) Mito genes are symbols here.
  visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "^MT-", assay = "Spatial")

  # QC plots BEFORE filtering
  vln <- VlnPlot(visium, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                 pt.size = 0.1, ncol = 3) & NoLegend()
  sp  <- SpatialFeaturePlot(visium, features = "nCount_Spatial") + theme(legend.position = "right")
  ggsave(file.path(out_dir, paste0("vlnplot_QC_sample", i, ".png")), vln, width = 10, height = 4)
  ggsave(file.path(out_dir, paste0("spatial_nCount_sample", i, ".png")), sp,  width = 8, height = 6)

  n_before <- ncol(visium)
  visium <- subset(visium, subset = nFeature_Spatial >= QC$min_features &
                                     nCount_Spatial   >= QC$min_counts &
                                     percent.mt        < QC$max_mt_percent)
  message("  QC: kept ", ncol(visium), " / ", n_before, " spots ",
          "(>=", QC$min_features, " genes, >=", QC$min_counts, " counts, <",
          QC$max_mt_percent, "% mito)")

  # ---- Normalize (SCTransform, regressing out mito) + reduction + clustering ----
  visium <- SCTransform(visium, assay = "Spatial", vars.to.regress = "percent.mt", verbose = FALSE)
  visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
  visium <- FindNeighbors(visium, reduction = "pca", dims = 1:p$dims)
  visium <- FindClusters(visium, resolution = p$resolution, verbose = FALSE)
  visium <- RunUMAP(visium, reduction = "pca", dims = 1:p$dims)

  # ---- Optional QC-driven cluster removal (e.g. empty cluster 8 in 1160920F) ----
  # Uses idents-based subsetting, which is robust to the factor levels of
  # seurat_clusters (the original `seurat_clusters != 8` comparison was fragile).
  if (length(p$remove_clusters) > 0) {
    Idents(visium) <- "seurat_clusters"
    keep <- setdiff(levels(visium), as.character(p$remove_clusters))
    visium <- subset(visium, idents = keep)
    message("  removed clusters: ", paste(p$remove_clusters, collapse = ", "))
  }

  # ---- Cluster visualization: UMAP + spatial ----
  clus <- DimPlot(visium, reduction = "umap", label = TRUE) +
          SpatialDimPlot(visium, label = TRUE, label.size = 3)
  ggsave(file.path(out_dir, paste0("clustering_plots_sample", i, ".png")), clus, width = 12, height = 6)

  # ---- Spatially variable features (Moran's I on top 1000 variable genes) ----
  # Restricted to 1000 genes because Moran's I over all genes is prohibitively slow.
  visium <- FindSpatiallyVariableFeatures(
    visium, assay = "SCT",
    features = VariableFeatures(visium)[1:1000],
    selection.method = "moransi"
  )
  top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)
  svf <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
  ggsave(file.path(out_dir, paste0("spatially_variable_features_s", i, ".png")), svf, width = 10, height = 6)

  # ---- Save processed object ----
  saveRDS(visium, processed_rds(id))
  message("  saved: ", processed_rds(id))
}

message("Step 1 complete.")
