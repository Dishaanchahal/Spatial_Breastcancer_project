# Part -2 
# Single cell script
setwd("/scratch/chahal.d/Spatial_transcriptomics/Wu_etal_2021_BRCA_scRNASeq")
library(Seurat)
library(ggplot2)

# Checking what we have after unzipping the files
meta_check <- read.csv("metadata.csv",row.names = 1)
head(meta_check)

# from the metadata there is a column called cell type major, checking what is there 
table(meta_check$celltype_major)

# Okay so it has all the major cell types - B-cells,CAFs, Cancer Epithelial, ENdothelial, Myeloid, Normal Epithelial, Plasmablasts, PVL, T-cells

# Building the reference seurat object
ref_counts <- Read10X(
  data.dir = ".",
  gene.column = 1
)

ref <- CreateSeuratObject(counts = ref_counts, meta.data = meta_check)
ref

# checking if QC was done
VlnPlot(ref, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
        pt.size = 0, ncol = 3)

# QC looks right
# Normalising using SCTransform the same way that we did for spatial data
# this will help with label transfer
ref <- SCTransform(ref, ncells = 3000, verbose = FALSE)
ref <- RunPCA(ref, verbose = FALSE)

ref <- RunUMAP(ref, dims = 1:30, verbose = FALSE)

p <- DimPlot(ref, group.by = "celltype_major", label = TRUE, repel = TRUE, raster=FALSE) + 
  NoLegend() +
  ggtitle("scRNA-seq Reference - Cell Types")
p

ggsave("reference_umap.png", plot = p, width = 10, height = 8, dpi = 300)

# save the seurat object
saveRDS(ref, "ref_processed.rds")

