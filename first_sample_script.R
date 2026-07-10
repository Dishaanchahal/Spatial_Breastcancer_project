# Trial script to learn Spatial Transcriptomics
# 17th April 2026
# Pipeline - Seurat's Spatial Transcriptomics

setwd("/scratch/chahal.d/Spatial_transcriptomics")
#Loading required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#############################Sample1
#############################

# Loading the datasets
# reading count matrix
counts <- Read10X("samples/1142243F/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/1142243F/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["1142243F"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

print(plot1)


plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2
# saving the plots
ggsave("vlnplot_nCount_sample1.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample1.png", plot = plot2, width = 8, height = 6)

# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)

# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)

# Find Clusters
visium <- FindClusters(visium, verbose = FALSE)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)
p1 + p2
png("clustering_plots_sample1.png", width = 1600, height = 800)
print(p1+p2)

dev.off()


#FInding spatially variable features
# Selecting top highly spatially variable genes
# Selecting high expression of genes clustering in a specific region, rathar than randomly
# This Moran's selection method checks that the gene X is surrounded by the Gene X 
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

# We plot top 6 spatially variable genes
plot<-SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
print(plot)

ggsave("spatial_features.png", plot = plot, width = 10, height = 6)
saveRDS(visium, file = "visium_sample1.rds")

############################# Sample2
############################# remember to run it again, cluster 8 removal failed

# Loading the datasets
# reading count matrix
counts <- Read10X("samples/1160920F/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/1160920F/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["1160920F"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

plot1
dev.off()
plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2

# saving the plots
ggsave("vlnplot_nCount_sample2.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample2.png", plot = plot2, width = 8, height = 6)



# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)

# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)

# Find Clusters
visium <- FindClusters(visium, verbose = FALSE)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)

plot <- p1 + p2
plot
ggsave("clustering_plots_sample2.png", plot = plot, width = 12, height = 6)

#Now in my UMAP cluster 8 and 10 seemed very off to me
# 8 Because it was clustered very away from all other clusters
# 10 because it was scattered
# Now I create a violin plot to check the number of cells and genes in those clusters
VlnPlot(visium, 
        features = c("nCount_Spatial", "nFeature_Spatial"), 
        pt.size = 0.1,
        idents = c("8", "10")
)

# Cluster 8 is clearly very much apart from all other and it has less numbers of genes and cells hence we should remove it, because
# Probably it is empty or nearly empty cells
# Cluster 10 actually has a lot of cells and genes so that is real, probably shares something with a lot of type of cells that is why the cluster is all over the place

#subsetting clusters and removing  8
visium <- subset(visium, subset = seurat_clusters != 8)

#FInding spatially variable features
# Selecting top highly spatially variable genes
# Selecting high expression of genes clustering in a specific region, rathar than randomly
# This Moran's selection method checks that the gene X is surrounded by the Gene X 
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

plot <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
plot
ggsave("spatially_variable_features_s2.png", plot = plot, width = 10, height = 6)
saveRDS(visium, file = "visium_sample2.rds")

######################## Sample 3
########################
# Loading the datasets
# reading count matrix
counts <- Read10X("samples/CID4290/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/CID4290/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["CID4290"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

plot1
dev.off()
plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2
dev.off()
# saving the plots
ggsave("vlnplot_nCount_sample3.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample3.png", plot = plot2, width = 8, height = 6)



# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)

# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:20)

# Find Clusters, # tried different resolutions, resolution 0.6 is giving the best clustering
visium <- FindClusters(visium, resolution = 0.6)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:20)
ElbowPlot(visium)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)

plot <- p1 + p2
plot

# Checking if CLuster 7 is real because it is scattered everywhere
Vlnplt<-VlnPlot(visium, features = c("nCount_Spatial", "nFeature_Spatial"),
        pt.size = 0.1, idents = c("7"))


ggsave("clustering_plots_sample3.png", plot = plot, width = 12, height = 6)
ggsave("cluster7_sample3_qc.png", plot = Vlnplt, width = 12, height =6)

#FInding spatially variable features
# Selecting top highly spatially variable genes
# Selecting high expression of genes clustering in a specific region, rathar than randomly
# This Moran's selection method checks that the gene X is surrounded by the Gene X 
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

plot <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
plot
ggsave("spatially_variable_features_s3.png", plot = plot, width = 10, height = 6)

saveRDS(visium, file = "visium_sample3.rds")

######################## Sample 4
########################
# Loading the datasets
# reading count matrix
counts <- Read10X("samples/CID4465/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/CID4465/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["CID4465"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

plot1

plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2

# saving the plots
ggsave("vlnplot_nCount_sample4.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample4.png", plot = plot2, width = 8, height = 6)



# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
ElbowPlot(visium)

# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:20)

# Find Clusters
visium <- FindClusters(visium, resolution = 0.5)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:20)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)

plot <- p1 + p2
plot
ggsave("clustering_plots_sample4.png", plot = plot, width = 12, height = 6)

#spatial
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

plot <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
plot
ggsave("spatially_variable_features_s4.png", plot = plot, width = 10, height = 6)
saveRDS(visium, file = "visium_sample4.rds")

#################### Sample 5
####################
# Loading the datasets
# reading count matrix
counts <- Read10X("samples/CID44971/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/CID44971/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["CID44971"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

plot1

plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2

# saving the plots
ggsave("vlnplot_nCount_sample5.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample5.png", plot = plot2, width = 8, height = 6)

# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)

ElbowPlot(visium)
# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:20)

# Find Clusters
visium <- FindClusters(visium, resolution = 0.5)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)

plot <- p1 + p2
plot
ggsave("clustering_plots_sample5.png", plot = plot, width = 12, height = 6)

#spatial
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

plot <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
plot
ggsave("spatially_variable_features_s5.png", plot = plot, width = 10, height = 6)
saveRDS(visium, file = "visium_sample5.rds")

#################### Sample 6
####################
# Loading the datasets
# reading count matrix
counts <- Read10X("samples/CID4535/filtered_feature_bc_matrix",gene.column=1)

visium <- CreateSeuratObject(counts, assay = "Spatial")
image <- Read10X_Image("samples/CID4535/spatial")
image <- image[Cells(visium)]
DefaultAssay(image) <- "Spatial"
visium[["CID4535"]] <- image

# sanity check
visium
head(visium@meta.data)

# QC Plots 
# This violin plot is to check the distibution of UMI per spot.
# X axis is sample
# Y axis is count per spot

# Second plot is no. of UMI's colored on the sample image. 
# Ideally TUmor and immune cells will have higher number of UMIs.
# Low count spots are fats
# If we see high counts in empty regions, that means something is wrong
plot1<-VlnPlot(visium, features ="nCount_Spatial",pt.size = 0.1) + NoLegend()

plot1

plot2<- SpatialFeaturePlot(visium, features ="nCount_Spatial")+ theme(legend.position="right")

plot2

# saving the plots
ggsave("vlnplot_nCount_sample6.png", plot = plot1, width = 6, height = 4)
ggsave("spatial_nCount_sample6.png", plot = plot2, width = 8, height = 6)

# Normalize with SCTransform 
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Running PCA
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
ElbowPlot(visium)
# Find neighbours
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:20)

# Find Clusters
visium <- FindClusters(visium, resolution = 0.5)

# Running UMAP
visium <- RunUMAP(visium, reduction = "pca", dims = 1:20)

# Now visualizing this
p1 <- DimPlot(visium, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(visium, label = TRUE, label.size = 3)

plot <- p1 + p2
plot
ggsave("clustering_plots_sample6.png", plot = plot, width = 12, height = 6)


# analysing cluster 0
VlnPlot(visium, features = c("nCount_Spatial", "nFeature_Spatial"),
        pt.size = 0.1, idents = c("1", "2"))


p <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))

ggsave("spatial_QC_Sample6.png", plot = p, width = 14, height = 8, dpi = 300)


### Also did find marker on CLuster 2, found it has tumor epithelial cells 
markers_c2 <- FindMarkers(visium, ident.1 = 2)
head(markers_c2, 20)

# It gave all tumor genes t

#spatial
visium <- FindSpatiallyVariableFeatures(
  visium, 
  assay = "SCT", 
  features = VariableFeatures(visium)[1:1000],
  selection.method = "moransi"
)

# Calculating Moransi selection on ALL genes will take forever that is why we are doing it on top 1000 genes
top.features <- head(SpatiallyVariableFeatures(visium, selection.method = "moransi"), 6)

plot <- SpatialFeaturePlot(visium, features = top.features, ncol = 3, alpha = c(0.1, 1))
plot
ggsave("spatially_variable_features_s6.png", plot = plot, width = 10, height = 6)

saveRDS(visium, file = "visium_sample6.rds")

