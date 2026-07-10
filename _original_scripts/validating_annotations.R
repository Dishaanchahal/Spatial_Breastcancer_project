# Validating pathologists annotations
library(pheatmap)
library(pheatmap)

# build validation heatmap for all samples as example
for(sample_name in c("s1","s2", "s3", "s4", "s5", "s6")) {
  obj <- get(sample_name)
  meta <- get(paste0(sample_name, "_meta"))
  df <- as.data.frame(t(GetAssayData(obj, assay = "predictions")))
  df$Classification <- meta[rownames(df), "Classification"]
  agg <- aggregate(. ~ Classification, data = df, FUN = mean)
  rownames(agg) <- agg$Classification
  agg$Classification <- NULL
  agg$max <- NULL
  pheatmap(agg,
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("navy", "white", "red"))(50),
           main = paste("Mean Prediction Score -", sample_name),
           filename = paste0("validation_heatmap_", sample_name, ".png"),
           width = 10, height = 6)
}
