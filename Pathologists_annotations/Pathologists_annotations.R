# pathologists annotations

# Before the label transfer and spatial transcriptomics, a pathologist looked at the H&E( stained images) and labelled
# spots as cancer, stromal,  lympocytes etc. 

# So we can compare spots which have high lymphocytes according to the pathologist, do they have high T/B cells?

library(Seurat)
library(ggplot2)
library(patchwork)

setwd("/scratch/chahal.d/Spatial_transcriptomics")

#Metadata
s1_meta <- read.csv("metadata/1142243F_metadata.csv", row.names = 1)
s2_meta <- read.csv("metadata/1160920F_metadata.csv", row.names = 1)
s3_meta <- read.csv("metadata/CID4290_metadata.csv", row.names = 1)
s4_meta <- read.csv("metadata/CID4465_metadata.csv", row.names = 1)
s5_meta <- read.csv("metadata/CID44971_metadata.csv", row.names = 1)
s6_meta <- read.csv("metadata/CID4535_metadata.csv", row.names = 1)

cat(" 1142243F (TNBC) ")
print(table(s1_meta$Classification))

cat(" 1160920F (TNBC) ")
print(table(s2_meta$Classification))

cat("CID4290 (ER+) ")
print(table(s3_meta$Classification))

cat(" CID4465 (TNBC) ")
print(table(s4_meta$Classification))

cat("CID44971 (TNBC)")
print(table(s5_meta$Classification))

cat(" CID4535 (ER+) ")
print(table(s6_meta$Classification))

# Loading the Spatial objects
s1 <- readRDS("sample1_1142243F_labeled.rds")
s2 <- readRDS("sample2_1160920F_labeled.rds")
s3 <- readRDS("sample3_CID4290_labeled.rds")
s4 <- readRDS("sample4_CID4465_labeled.rds")
s5 <- readRDS("sample5_CID44971_labeled.rds")
s6 <- readRDS("sample6_CID4535_labeled.rds")

# Adding the pathologists annotations
s1 <- AddMetaData(s1, s1_meta)
s2 <- AddMetaData(s2, s2_meta)
s3 <- AddMetaData(s3, s3_meta)
s4 <- AddMetaData(s4, s4_meta)
s5 <- AddMetaData(s5, s5_meta)
s6 <- AddMetaData(s6, s6_meta)

# Plot pathologists annotations spatially 
p1 <- SpatialDimPlot(s1, group.by = "Classification") + ggtitle("1142243F TNBC")
p2 <- SpatialDimPlot(s2, group.by = "Classification") + ggtitle("1160920F TNBC")
p3 <- SpatialDimPlot(s3, group.by = "Classification") + ggtitle("CID4290 ER+")
p4 <- SpatialDimPlot(s4, group.by = "Classification") + ggtitle("CID4465 TNBC")
p5 <- SpatialDimPlot(s5, group.by = "Classification") + ggtitle("CID44971 TNBC")
p6 <- SpatialDimPlot(s6, group.by = "Classification") + ggtitle("CID4535 ER+")

p1
p2
p3
p4
p5
p6

ggsave("pathologist_annotations_all.png", 
       plot = p1 + p2 + p3 + p4 + p5 + p6,
       width = 18, height = 12, dpi = 300)

# Comparing these to what we found, the pathologists annotations align well. Our predictions from spatial data are validated