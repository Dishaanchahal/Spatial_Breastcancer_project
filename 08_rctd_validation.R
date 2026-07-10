# Option 3 (validate): independent deconvolution with RCTD (spacexr),
# to cross-check the cell2location composition estimates.
suppressPackageStartupMessages({ library(Seurat); library(spacexr); library(Matrix) })
source("00_config.R")
OUT <- file.path(PROJECT_ROOT, "Validation_RCTD"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(0)

ref <- readRDS(REF_RDS)
# subsample per major cell type (match cell2location: <=2500/type) for tractability
keep <- unlist(lapply(split(colnames(ref), ref$celltype_major),
                      function(cc) if (length(cc) > 2500) sample(cc, 2500) else cc))
ref <- subset(ref, cells = keep)
ref_counts <- GetAssayData(ref, assay = "RNA", slot = "counts")
cell_types <- as.factor(ref$celltype_major); names(cell_types) <- colnames(ref)
nUMI_ref   <- colSums(ref_counts);           names(nUMI_ref)   <- colnames(ref)
reference  <- Reference(ref_counts, cell_types, nUMI_ref)

comp <- list()
for (id in SAMPLES) {
  message("=== RCTD | ", id, " ===")
  obj    <- readRDS(processed_rds(id))
  counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
  coords <- as.data.frame(GetTissueCoordinates(obj))
  coords <- coords[, sapply(coords, is.numeric)][, 1:2]; colnames(coords) <- c("x", "y")
  puck   <- SpatialRNA(coords, counts, colSums(counts))
  rc     <- create.RCTD(puck, reference, max_cores = 4)
  rc     <- run.RCTD(rc, doublet_mode = "full")
  w      <- as.matrix(rc@results$weights)
  w      <- sweep(w, 1, rowSums(w), "/")
  write.csv(w, file.path(OUT, paste0("rctd_weights_", id, ".csv")))
  comp[[id]] <- colMeans(w)
}
df <- as.data.frame(do.call(rbind, comp)); df$sample <- rownames(df); df$subtype <- SUBTYPES[df$sample]
write.csv(df, file.path(OUT, "rctd_composition.csv"), row.names = FALSE)
message("RCTD_DONE")
