# Compare cell2location vs RCTD per-sample composition. Base R only (robust to
# missing ggplot2/reshape2), cairo PNG device (headless-safe).
source("00_config.R")
options(bitmapType = "cairo")
OUT  <- file.path(PROJECT_ROOT, "Validation_RCTD")
c2l  <- read.csv(file.path(PROJECT_ROOT, "cell2location", "results", "tme_composition_c2l.csv"), check.names = FALSE)
rctd <- read.csv(file.path(OUT, "rctd_composition.csv"), check.names = FALSE)
rownames(c2l) <- c2l$sample; rownames(rctd) <- rctd$sample
cts <- CELL_TYPES
samples <- intersect(rownames(c2l), rownames(rctd))

m <- do.call(rbind, lapply(samples, function(s) do.call(rbind, lapply(cts, function(ct)
      data.frame(sample = s, cell_type = ct, c2l = c2l[s, ct], rctd = rctd[s, ct])))))
r <- cor(m$c2l, m$rctd)
write.csv(m, file.path(OUT, "c2l_vs_rctd_merged.csv"), row.names = FALSE)

cols <- setNames(rainbow(length(cts)), cts)
lim  <- c(0, max(m$c2l, m$rctd))
png(file.path(OUT, "c2l_vs_rctd_scatter.png"), width = 1700, height = 1500, res = 200, type = "cairo")
plot(m$c2l, m$rctd, col = cols[m$cell_type], pch = 19, cex = 1.4, xlim = lim, ylim = lim,
     xlab = "cell2location proportion", ylab = "RCTD proportion",
     main = sprintf("cell2location vs RCTD composition (Pearson r = %.2f)", r))
abline(0, 1, lty = 2, col = "grey50")
legend("topleft", legend = cts, col = cols, pch = 19, cex = 0.8, bty = "n")
invisible(dev.off())
cat(sprintf("Pearson r (cell2location vs RCTD) = %.3f\n", r))
message("COMPARE_DONE")
