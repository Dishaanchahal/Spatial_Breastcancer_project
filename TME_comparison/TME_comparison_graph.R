# TME comparison graph 

library(reshape2)
library(ggplot2)

# build summary dataframe
tme_summary <- data.frame(
  Sample = c("1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"),
  Cancer_Epithelial = c(0.318, 0.611, 0.934, 0.168, 0.319, 0.425),
  T_cells           = c(0.143, 0.027, 0.000, 0.000, 0.032, 0.312),
  CAFs              = c(0.194, 0.008, 0.062, 0.514, 0.250, 0.064),
  Normal_Epithelial = c(0.179, 0.062, 0.000, 0.016, 0.220, 0.000),
  Plasmablasts      = c(0.101, 0.219, 0.000, 0.245, 0.099, 0.134),
  Myeloid           = c(0.023, 0.024, 0.000, 0.031, 0.019, 0.044),
  B_cells           = c(0.024, 0.048, 0.000, 0.007, 0.015, 0.014),
  Endothelial       = c(0.009, 0.000, 0.001, 0.007, 0.014, 0.003)
)

# melt to long format
tme_long <- melt(tme_summary, id.vars = "Sample", 
                 variable.name = "Cell_Type", 
                 value.name = "Score")

# plot

p <- ggplot(tme_long, aes(x = Sample, y = Score, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TME Composition Across 6 Breast Cancer Samples",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent)
p

ggsave("tme_comparison.png", plot = p, width = 10, height = 6, dpi = 300)

# subtypes added
# update tme_summary with subtype
tme_summary$Subtype <- c("TNBC", "TNBC", "ER+", "TNBC", "TNBC", "ER+")

# melt again
tme_long <- melt(tme_summary, id.vars = c("Sample", "Subtype"),
                 variable.name = "Cell_Type",
                 value.name = "Score")

# plot with subtype on x axis
p <- ggplot(tme_long, aes(x = Sample, y = Score, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Subtype, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(title = "TME Composition by Breast Cancer Subtype",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent)

ggsave("tme_comparison_subtype.png", plot = p, width = 12, height = 6, dpi = 300)