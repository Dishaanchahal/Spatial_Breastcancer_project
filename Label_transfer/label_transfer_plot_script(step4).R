# VIsualising the label transfers

setwd("/scratch/chahal.d/Spatial_transcriptomics")
library(Seurat)
library(ggplot2)

##### sample 1
s1 <- readRDS("sample1_1142243F_labeled.rds")

# set predictions as active assay
DefaultAssay(s1) <- "predictions"

# plot all 9 cell types spatially
p <- SpatialFeaturePlot(s1, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)

p
ggsave("label_transfer_s1.png", plot = p, width = 14, height = 14, dpi = 300)


# sample 2
s2 <- readRDS("sample2_1160920F_labeled.rds")
DefaultAssay(s2) <- "predictions"

# check scores
rowMeans(GetAssayData(s2))

# plot
p <- SpatialFeaturePlot(s2, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)

ggsave("label_transfer_s2.png", plot = p, width = 14, height = 14, dpi = 300)




########## sample 3
s3 <- readRDS("sample3_CID4290_labeled.rds")
DefaultAssay(s3) <- "predictions"

rowMeans(GetAssayData(s3))

p <- SpatialFeaturePlot(s3, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)
p

ggsave("label_transfer_s3.png", plot = p, width = 14, height = 14, dpi = 300)



############ sample 4
s4 <- readRDS("sample4_CID4465_labeled.rds")
DefaultAssay(s4) <- "predictions"

rowMeans(GetAssayData(s4))

p <- SpatialFeaturePlot(s4, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)
p
ggsave("label_transfer_s4.png", plot = p, width = 14, height = 14, dpi = 300)

############### sample5
s5 <- readRDS("sample5_CID44971_labeled.rds")
DefaultAssay(s5) <- "predictions"

rowMeans(GetAssayData(s5))

p <- SpatialFeaturePlot(s5, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)
p
ggsave("label_transfer_s5.png", plot = p, width = 14, height = 14, dpi = 300)

#################### sample6

s6 <- readRDS("sample6_CID4535_labeled.rds")
DefaultAssay(s6) <- "predictions"

rowMeans(GetAssayData(s6))

p <- SpatialFeaturePlot(s6, 
                        features = c("Cancer Epithelial", "T-cells", "Myeloid", 
                                     "CAFs", "Endothelial", "B-cells", 
                                     "Plasmablasts", "PVL", "Normal Epithelial"),
                        ncol = 3, 
                        alpha = c(0.1, 1),
                        min.cutoff = 0,
                        max.cutoff = 0.5
) & scale_fill_gradientn(
  colors = c("navy", "blue", "cyan", "yellow", "red"),
  limits = c(0, 0.5)
)
p
ggsave("label_transfer_s6.png", plot = p, width = 14, height = 14, dpi = 300)
