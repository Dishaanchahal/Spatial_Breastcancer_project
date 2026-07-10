# Option 3 (validate): independent deconvolution with RCTD (spacexr).
# Reads raw counts/coords straight from 10x flat files (no Seurat dependency).
suppressPackageStartupMessages({ library(spacexr); library(Matrix) })
source("00_config.R")
OUT <- file.path(PROJECT_ROOT, "Validation_RCTD"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(0)

is_gz <- function(path) { con <- file(path, "rb"); m <- readBin(con, "raw", 2); close(con); identical(as.integer(m), c(31L, 139L)) }
gunzip_tmp <- function(path) {
  tmp <- tempfile(); inp <- gzfile(path, "rb"); out <- file(tmp, "wb")
  repeat { b <- readBin(inp, "raw", 1e7); if (length(b) == 0) break; writeBin(b, out) }
  close(inp); close(out); tmp
}
read_lines2 <- function(path) { p <- if (is_gz(path)) gunzip_tmp(path) else path; readLines(p) }
read_mtx2   <- function(path) { p <- if (is_gz(path)) gunzip_tmp(path) else path; as(Matrix::readMM(p), "CsparseMatrix") }
load_counts <- function(dir) {
  bc <- read_lines2(file.path(dir, "barcodes.tsv.gz"))
  ft <- read_lines2(file.path(dir, "features.tsv.gz")); genes <- sub("\t.*", "", ft)
  m  <- read_mtx2(file.path(dir, "matrix.mtx.gz"))
  if (nrow(m) == length(bc) && ncol(m) == length(genes)) m <- t(m)
  rownames(m) <- make.unique(genes); colnames(m) <- bc
  round(m)
}

ref_counts <- load_counts(DIR_REF)
meta   <- read.csv(file.path(DIR_REF, "metadata.csv"), row.names = 1)
common <- intersect(colnames(ref_counts), rownames(meta))
ref_counts <- ref_counts[, common]; ct <- meta[common, "celltype_major"]
idx <- unlist(lapply(split(seq_along(ct), ct), function(ii) if (length(ii) > 2500) sample(ii, 2500) else ii))
ref_counts <- ref_counts[, idx]; ct <- ct[idx]
cell_types <- as.factor(ct); names(cell_types) <- colnames(ref_counts)
nUMI_ref <- colSums(ref_counts); names(nUMI_ref) <- colnames(ref_counts)
reference <- Reference(ref_counts, cell_types, nUMI_ref)
message("reference: ", ncol(ref_counts), " cells x ", nrow(ref_counts), " genes")

comp <- list()
for (id in SAMPLES) {
  message("=== RCTD | ", id, " ===")
  counts <- load_counts(file.path(DIR_SAMPLES, id, "filtered_feature_bc_matrix"))
  pos <- read.csv(file.path(DIR_SAMPLES, id, "spatial", "tissue_positions_list.csv"), header = FALSE)
  rownames(pos) <- pos[, 1]
  coords <- pos[colnames(counts), c(6, 5)]; colnames(coords) <- c("x", "y")
  ok <- complete.cases(coords); counts <- counts[, ok]; coords <- coords[ok, ]
  nUMI <- colSums(counts); keep <- nUMI > 0
  counts <- counts[, keep]; coords <- coords[keep, ]; nUMI <- nUMI[keep]
  puck <- SpatialRNA(coords, counts, nUMI)
  rc <- create.RCTD(puck, reference, max_cores = 4)
  rc <- run.RCTD(rc, doublet_mode = "full")
  w  <- as.matrix(rc@results$weights); w <- sweep(w, 1, rowSums(w), "/")
  write.csv(w, file.path(OUT, paste0("rctd_weights_", id, ".csv")))
  comp[[id]] <- colMeans(w)
}
df <- as.data.frame(do.call(rbind, comp)); df$sample <- rownames(df); df$subtype <- SUBTYPES[df$sample]
write.csv(df, file.path(OUT, "rctd_composition.csv"), row.names = FALSE)
message("RCTD_DONE")
