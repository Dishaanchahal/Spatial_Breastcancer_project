# Compare cell2location vs RCTD per-sample composition (concordance).
suppressPackageStartupMessages({ library(ggplot2); library(reshape2) })
source("00_config.R")
OUT  <- file.path(PROJECT_ROOT, "Validation_RCTD")
c2l  <- read.csv(file.path(PROJECT_ROOT, "cell2location", "results", "tme_composition_c2l.csv"), check.names = FALSE)
rctd <- read.csv(file.path(OUT, "rctd_composition.csv"), check.names = FALSE)
cts  <- CELL_TYPES
mlt  <- function(d) reshape2::melt(d[, c("sample", intersect(cts, colnames(d)))],
                                   id.vars = "sample", variable.name = "cell_type", value.name = "prop")
wide <- merge(mlt(c2l), mlt(rctd), by = c("sample", "cell_type"), suffixes = c("_c2l", "_rctd"))
r <- cor(wide$prop_c2l, wide$prop_rctd)
p <- ggplot(wide, aes(prop_c2l, prop_rctd, color = cell_type)) +
  geom_abline(linetype = "dashed", color = "grey60") + geom_point(size = 2.5) +
  labs(x = "cell2location proportion", y = "RCTD proportion",
       title = sprintf("cell2location vs RCTD composition (Pearson r = %.2f)", r)) +
  theme_classic()
ggsave(file.path(OUT, "c2l_vs_rctd_scatter.png"), p, width = 8, height = 6, dpi = 150)
write.csv(wide, file.path(OUT, "c2l_vs_rctd_merged.csv"), row.names = FALSE)
cat(sprintf("Pearson r (cell2location vs RCTD) = %.3f\n", r))
message("COMPARE_DONE")
