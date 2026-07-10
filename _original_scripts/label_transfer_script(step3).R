# Label transfer

library(Seurat)
library(ggplot2)
library(patchwork)

setwd("/scratch/chahal.d/Spatial_transcriptomics")

# read the Single cell  rds to transfer labels from
ref <- readRDS("Wu_etal_2021_BRCA_scRNASeq/ref_processed.rds")

# read the Spatial rds to transfer the labels to 
s1 <- readRDS("Processed_objects_spatial/sample_1/visium_sample1.rds")
s2 <- readRDS("Processed_objects_spatial/sample_2/visium_sample2.rds")
s3 <- readRDS("Processed_objects_spatial/sample_3/visium_sample3.rds")
s4 <- readRDS("Processed_objects_spatial/sample_4/visium_sample4.rds")
s5 <- readRDS("Processed_objects_spatial/sample_5/visium_sample5.rds")
s6 <- readRDS("Processed_objects_spatial/sample_6/visium_sample6.rds")

# Function to transfer the lables
run_label_transfer <- function(spatial_obj, sample_name) {
  
  cat("Processing", sample_name, "\n")
  
  anchors <- FindTransferAnchors(
    reference = ref,
    query = spatial_obj,
    normalization.method = "SCT",
    verbose = FALSE
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref$celltype_major,
    prediction.assay = TRUE,
    weight.reduction = spatial_obj[["pca"]],
    dims = 1:30
  )
  
  spatial_obj[["predictions"]] <- predictions
  return(spatial_obj)
}

# run label transfer on all samples

s1 <- run_label_transfer(s1, "1142243F")
s2 <- run_label_transfer(s2, "1160920F")
s3 <- run_label_transfer(s3, "CID4290")
s4 <- run_label_transfer(s4, "CID4465")
s5 <- run_label_transfer(s5, "CID44971")
s6 <- run_label_transfer(s6, "CID4535")

# save label transferred rds
saveRDS(s1, "sample1_1142243F_labeled.rds")
saveRDS(s2, "sample2_1160920F_labeled.rds")
saveRDS(s3, "sample3_CID4290_labeled.rds")
saveRDS(s4, "sample4_CID4465_labeled.rds")
saveRDS(s5, "sample5_CID44971_labeled.rds")
saveRDS(s6, "sample6_CID4535_labeled.rds")
