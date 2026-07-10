# ============================================================================
# 07_tme_comparison.R  —  Tumor-microenvironment composition across samples
#
# Refactor of TME_comparison/TME_comparison_graph.R. The original hardcoded the
# per-sample mean prediction scores. Here they are COMPUTED from the labeled
# objects, so the figure always reflects the current analysis and is fully
# reproducible from the data.
# ============================================================================

source("00_config.R")
suppressPackageStartupMessages({ library(Seurat); library(reshape2); library(ggplot2) })

# ---- Compute mean prediction score per cell type, per sample ----
score_list <- lapply(SAMPLES, function(id) {
  obj    <- readRDS(labeled_rds(id))
  scores <- rowMeans(GetAssayData(obj, assay = "predictions"))
  scores <- scores[setdiff(names(scores), "max")]  # drop TransferData 'max' row
  scores
})
names(score_list) <- SAMPLES

# samples x cell types
mat <- do.call(rbind, lapply(score_list, function(s) s[CELL_TYPES]))
colnames(mat) <- CELL_TYPES

tme_summary <- data.frame(Sample = SAMPLES, mat, check.names = FALSE)
tme_summary$Subtype <- SUBTYPES[SAMPLES]
write.csv(tme_summary, file.path(DIR_TME, "tme_summary_computed.csv"), row.names = FALSE)

# ---- Long format ----
tme_long <- melt(tme_summary, id.vars = c("Sample", "Subtype"),
                 variable.name = "Cell_Type", value.name = "Score")

# ---- Plot 1: composition across all samples ----
p1 <- ggplot(tme_long, aes(x = Sample, y = Score, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TME Composition Across 6 Breast Cancer Samples",
       subtitle = "Label-transfer prediction scores (approximate; see 03b_deconvolution_rctd.R)",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent)
ggsave(file.path(DIR_TME, "tme_comparison.png"), plot = p1, width = 10, height = 6, dpi = 300)

# ---- Plot 2: faceted by subtype (TNBC vs ER+) ----
p2 <- ggplot(tme_long, aes(x = Sample, y = Score, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Subtype, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 12, face = "bold")) +
  labs(title = "TME Composition by Breast Cancer Subtype",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent)
ggsave(file.path(DIR_TME, "tme_comparison_subtype.png"), plot = p2, width = 12, height = 6, dpi = 300)

message("Step 7 complete. Wrote tme_summary_computed.csv + 2 figures.")
