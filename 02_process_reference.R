# ============================================================================
# 02_process_reference.R  —  Step 2: build the scRNA-seq reference
#
# Refactor of Wu_etal_2021_BRCA_scRNASeq/single_cell_processing(step2).R.
# Builds an SCT-normalized Seurat reference from the Wu et al. 2021 BRCA
# atlas so its cell-type labels can be transferred onto the Visium spots.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

meta <- read.csv(file.path(DIR_REF, "metadata.csv"), row.names = 1)
message("Reference cell types (celltype_major):")
print(table(meta$celltype_major))

ref_counts <- Read10X(data.dir = DIR_REF, gene.column = 1)
ref <- CreateSeuratObject(counts = ref_counts, meta.data = meta)

# QC sanity check (Wu et al. data is already filtered)
VlnPlot(ref, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0, ncol = 3)

# SCTransform to match the spatial normalization (required for label transfer)
ref <- SCTransform(ref, ncells = 3000, verbose = FALSE)
ref <- RunPCA(ref, verbose = FALSE)
ref <- RunUMAP(ref, dims = 1:30, verbose = FALSE)

p <- DimPlot(ref, group.by = "celltype_major", label = TRUE, repel = TRUE, raster = FALSE) +
  NoLegend() + ggtitle("scRNA-seq Reference - Cell Types")
ggsave(file.path(DIR_REF, "reference_umap.png"), plot = p, width = 10, height = 8, dpi = 300)

saveRDS(ref, REF_RDS)
message("Step 2 complete. Saved reference: ", REF_RDS)
