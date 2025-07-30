#!/usr/bin/env python3
# MIT License
#
# Copyright 2024 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import logging
import shutil
import os
from scipy import io
import pandas as pd
import anndata as ad
from tqdm import tqdm
import numpy as np


logger = logging.getLogger(__name__)
has_pigz_binary = shutil.which("pigz") is not None

try:
    from pygz import PigzFile
    has_pygz = True

except ImportError:
    PigzFile = None  # type: ignore
    has_pygz = False
    import gzip
    logger.warn("pygz not installed.  Using gzip instead, which may be slower.")

use_pygz = has_pygz and has_pigz_binary

def write_10x_mtx(adata: ad.AnnData, output_dir: str) -> None:
    """
    Write an AnnData object to 10x-style MTX format (compressed).

    Parameters
    ----------
    adata : AnnData
        The filtered AnnData object.
    output_dir : str
        Directory to write `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`.
    """
    logger.info(f"Writing 10x-format MTX files to: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # Write matrix.mtx.gz
    matrix_path = os.path.join(output_dir, "matrix.mtx.gz")
    io.mmwrite(matrix_path[:-3], adata.X.tocoo())
    os.rename(matrix_path[:-3], matrix_path)

    # Write barcodes.tsv.gz
    barcodes_path = os.path.join(output_dir, "barcodes.tsv.gz")
    adata.obs_names.to_series().to_csv(barcodes_path, index=False, header=False, sep="\t", compression="gzip")

    # Write features.tsv.gz
    features_path = os.path.join(output_dir, "features.tsv.gz")
    if "gene_id" in adata.var.columns and "gene_name" in adata.var.columns:
        features_df = adata.var[["gene_id", "gene_name"]].copy()
    else:
        features_df = pd.DataFrame({
            "gene_id": adata.var_names,
            "gene_name": adata.var_names,
        })

    features_df["feature_type"] = "Gene Expression"
    features_df.to_csv(features_path, index=False, header=False, sep="\t", compression="gzip")


def compress_matrix_file(matrix_txt_path: str) -> str:
    """
    Compress a text-based matrix file using pygz if available, otherwise gzip.

    Returns the path to the compressed file.
    """
    logger.info("Starting compression of matrix.mtx")
    matrix_gz_path = matrix_txt_path + ".gz"

    if use_pygz:
        logger.info("Using pygz for compression")
        with open(matrix_txt_path, "rt") as f_in, PigzFile(matrix_gz_path, "wt") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        logger.info("Using gzip for compression")
        with open(matrix_txt_path, "rt") as f_in, gzip.open(matrix_gz_path, "wt") as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(matrix_txt_path)
    logger.info("Ended compression of matrix.mtx")
    return matrix_gz_path

def write_10x_matrix_market_file_faster(
    matrix_path: str,
    zarr_chunk_paths: list[str],
    n_cells: int,
    n_genes: int,
    nnz: int,
    flush_every: int = 1_000_000
) -> None:
    """
    Write matrix.mtx.gz in 10x's Matrix Market format from Zarr chunks using buffered streaming.
    Writes to uncompressed .mtx first, then compresses to .gz.

    MTX format convention:
    - Rows = genes
    - Columns = cells
    - Indices are 1-based
    """
    uncompressed_path = matrix_path[:-3] if matrix_path.endswith(".gz") else matrix_path

    with open(uncompressed_path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{n_genes} {n_cells} {nnz}\n")  # Genes = rows, Cells = columns

        cell_offset = 0

        for path in tqdm(zarr_chunk_paths, desc="Pass 2: Writing matrix.mtx"):
            adata = ad.read_zarr(path)
            X = adata.X.tocoo()

            n_entries = X.data.shape[0]

            for i in range(0, n_entries, flush_every):
                # Correct interpretation of COO: row = cell, col = gene
                rows = X.col[i:i + flush_every] + 1                      # genes = rows in MTX
                cols = X.row[i:i + flush_every] + cell_offset + 1        # cells = columns in MTX
                vals = X.data[i:i + flush_every]

                row_strs = rows.astype(str)
                col_strs = cols.astype(str)
                val_strs = vals.astype(str)

                lines = np.char.add(np.char.add(row_strs, " "), np.char.add(col_strs, " "))
                lines = np.char.add(lines, val_strs)
                lines = np.char.add(lines, "\n")

                f.writelines(lines.tolist())

            cell_offset += adata.n_obs  # accumulate # of cells in each chunk

    # Gzip compress after writing
    compress_matrix_file(uncompressed_path)

def write_standard_matrix_market_file(
    matrix_path: str,
    zarr_chunk_paths: list[str],
    n_cells: int,
    n_genes: int,
    nnz: int,
    flush_every: int = 1_000_000
) -> None:
    """
    Write matrix.mtx.gz in standard Matrix Market format from Zarr chunks using buffered streaming.
    This version writes cells as rows and genes as columns, matching the AnnData .X layout.

    Writes to uncompressed .mtx first, then compresses to .gz.

    MTX format convention (standard):
    - Rows = cells
    - Columns = genes
    - Indices are 1-based
    """
    uncompressed_path = matrix_path[:-3] if matrix_path.endswith(".gz") else matrix_path

    with open(uncompressed_path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{n_cells} {n_genes} {nnz}\n")  # Cells = rows, Genes = columns

        row_offset = 0

        for path in tqdm(zarr_chunk_paths, desc="Pass 2: Writing matrix.mtx"):
            adata = ad.read_zarr(path)
            X = adata.X.tocoo()

            n_entries = X.data.shape[0]

            for i in range(0, n_entries, flush_every):
                rows = X.row[i:i + flush_every] + row_offset + 1  # cells = rows
                cols = X.col[i:i + flush_every] + 1                # genes = columns
                vals = X.data[i:i + flush_every]

                row_strs = rows.astype(str)
                col_strs = cols.astype(str)
                val_strs = vals.astype(str)

                lines = np.char.add(np.char.add(row_strs, " "), np.char.add(col_strs, " "))
                lines = np.char.add(lines, val_strs)
                lines = np.char.add(lines, "\n")

                f.writelines(lines.tolist())

            row_offset += adata.n_obs  # accumulate # of cells

    # Gzip compress after writing
    compress_matrix_file(uncompressed_path)
