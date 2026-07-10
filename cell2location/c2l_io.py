"""Robust 10x loader for this project.

Quirks handled:
- Some matrix/barcode/feature files are plain text despite a .gz extension
  (spatial samples), while others are truly gzipped (reference). We sniff the
  gzip magic bytes and open accordingly.
- features.tsv(.gz) here is a SINGLE column of gene symbols (not the standard
  2-3 column 10x format), so scanpy's read_10x_mtx cannot be used directly.
"""
import gzip
import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


def _open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def _read_lines(path):
    with _open_text(path) as fh:
        return [ln.rstrip("\n") for ln in fh]


def _read_mtx(path):
    with _open_text(path) as fh:
        m = scipy.io.mmread(fh)
    return sp.csr_matrix(m)


def load_10x(dir_path, feature_col=0):
    """Return AnnData (cells x genes) with raw counts."""
    barcodes = _read_lines(os.path.join(dir_path, "barcodes.tsv.gz"))
    feats = _read_lines(os.path.join(dir_path, "features.tsv.gz"))
    genes = [f.split("\t")[feature_col] for f in feats]

    M = _read_mtx(os.path.join(dir_path, "matrix.mtx.gz"))  # 10x = genes x cells
    if M.shape[0] == len(genes) and M.shape[1] == len(barcodes):
        X = M.T.tocsr()                    # -> cells x genes
    elif M.shape[0] == len(barcodes) and M.shape[1] == len(genes):
        X = M.tocsr()
    else:
        raise ValueError(f"matrix {M.shape} does not match "
                         f"{len(barcodes)} barcodes / {len(genes)} genes in {dir_path}")

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    return adata


def load_spatial_coords(spatial_dir):
    """tissue_positions_list.csv -> DataFrame indexed by barcode with x,y."""
    p = os.path.join(spatial_dir, "tissue_positions_list.csv")
    df = pd.read_csv(p, header=None)
    df.columns = ["barcode", "in_tissue", "array_row", "array_col",
                  "pxl_row", "pxl_col"][: df.shape[1]]
    df = df.set_index("barcode")
    return df
