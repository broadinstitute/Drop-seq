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
import os
import argparse
import shutil
import gc
import glob
import warnings

import anndata as ad
import numpy as np
import sys

from scipy.sparse import csr_matrix, issparse
from tqdm import tqdm
from scipy import io
import pandas as pd
from typing import List, Optional
import psutil
import tempfile
from scipy.sparse import hstack
import scipy.sparse as sp
import time

import dropseq.hdf5.filters as filters

import logging

from dropseq.hdf5.downsample_adata import downsample_anndata, compute_downsampling_counts
from dropseq.hdf5.mtx_writer import write_10x_matrix_market_file_faster, write_standard_matrix_market_file

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

# suppress the pandas FutureWarning about "empty or all-NA entries"
# We explicitly "sanitize" the DataFrame before writing to Zarr, so this is safe.
warnings.filterwarnings(
    "ignore",
    message=".*empty or all-NA entries is deprecated.*",
    category=FutureWarning,
)

def extract_raw_counts(adata: ad.AnnData, cast_to_int: bool = False) -> ad.AnnData:
    """
    Extracts the raw counts from an AnnData object, with an option to cast the data to integer type.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object containing the data.
    cast_to_int : bool, optional
        If True, casts the expression data to integer type. Default is False.

    Returns
    -------
    ad.AnnData
        A new AnnData object containing the raw counts, optionally cast to integers.

    Raises
    ------
    ValueError
        If the AnnData object does not have a `.raw` attribute.
    """
    if adata.raw is None:
        raise ValueError("The AnnData object does not have a `.raw` attribute.")

    # Extract the raw data
    raw_adata = adata.raw.to_adata()

    if cast_to_int:
        # Function to safely cast data to integers
        def safe_cast_to_int(X):
            if sp.issparse(X):
                return X.astype(np.int32)
            else:
                return np.rint(X).astype(np.int32)

        raw_adata.X = safe_cast_to_int(raw_adata.X)

    return raw_adata


def safe_log(msg: str) -> None:
    """
    Write a message using tqdm.write if tqdm is active, else use logging.info.
    This avoids triggering exceptions if tqdm.write is not safe.
    """
    if len(tqdm._instances) > 0:
        tqdm.write(msg)
    else:
        logger.info(msg)

def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and sanitize a DataFrame (e.g., `.obs` or `.var`) so it's safe for Zarr.

    - Converts columns with only numeric + "NA" string values to float
    - Fills NaNs with "NA" for object columns
    - Coerces object columns to string if needed

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to sanitize (e.g., `.obs` or `.var`).

    Returns
    -------
    pd.DataFrame
        A sanitized copy of the input DataFrame.
    """
    df = df.copy()

    for col in df.columns:
        series = df[col]

        # Special handling for categorical columns
        # If it’s a Categorical, cast to string and fill missing
        if isinstance(series.dtype, pd.CategoricalDtype):
            df[col] = series.astype(str).fillna("NA")
            continue

        if series.dtype == object:
            try:
                # Attempt numeric coercion, treating "NA" as missing
                converted = pd.to_numeric(series.replace("NA", np.nan))
                df[col] = converted
                continue
            except Exception:
                pass

        # TODO: this may be redundant.
        # Might be able to replace this and the special column handling above with:
        # df[col] = series.astype(str).fillna("NA")
        if series.isna().any():
            df[col] = series.fillna("NA").astype(str)
        elif series.dtype == object:
            df[col] = series.astype(str)

    return df

def log_memory(note: str = ""):
    gc.collect()
    mem_gb = psutil.Process(os.getpid()).memory_info().rss / 1024**3
    msg = f"Memory usage{f' ({note})' if note else ''}: {mem_gb:.2f} GB"
    safe_log(msg)

def as_dataframe(adata):
    df = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
    )
    return df

# def get_union_var_names(h5ad_files: List[str],) -> List[str]:
#     """
#     Loop through all .h5ad files to collect the full set of variable (gene) names.
#     Uses backed mode to avoid loading full data into memory.
#     Ensures all files are properly closed even if an error occurs.
#     """
#     var_names_union: Set[str] = set()
#
#     for fname in tqdm(h5ad_files, desc="Finding all gene symbols"):
#         adata = None
#         try:
#             adata = ad.read_h5ad(fname, backed='r')
#             var_names_union.update(adata.var_names)
#         except Exception as e:
#             tqdm.write(f"Failed to read {fname} for var names: {e}")
#         finally:
#             if adata is not None:
#                 adata.file.close()
#
#     log_memory("Finished gathering all gene names")
#     return sorted(var_names_union)

def get_union_var_df(h5ad_files: List[str],) -> pd.DataFrame:
    """
    Loop through all .h5ad files and collect the union of .var rows across files.
    Ensures that missing genes can be filled with appropriate metadata,
    and preserves index name if available.
    """
    var_dfs = []

    for fname in tqdm(h5ad_files, desc="Finding all gene metadata"):
        adata = None
        try:
            adata = ad.read_h5ad(fname, backed='r')
            df = adata.var.copy()
            if df.index.name is None:
                df.index.name = "gene_id"
            var_dfs.append(df)
        except Exception as e:
            tqdm.write(f"Failed to read {fname} for .var: {e}")
        finally:
            if adata is not None:
                adata.file.close()

    # Concatenate and drop duplicates based on index
    combined_var = pd.concat(var_dfs)

    if combined_var.index.name is None:
        combined_var.index.name = "gene_id"
    combined_var = combined_var[~combined_var.index.duplicated(keep="first")]
    combined_var = combined_var.sort_index()

    log_memory("Finished gathering all gene metadata")
    return combined_var

def measure_peak_memory(func, *args, **kwargs):
    """
    Measure peak memory before and after a function call using psutil, and log the result.
    Performs garbage collection before and after the function to improve accuracy.
    Returns the result of the function so it can be used downstream.
    """
    process = psutil.Process(os.getpid())

    gc.collect()
    mem_before = process.memory_info().rss

    start = time.time()
    result = func(*args, **kwargs)
    gc.collect()
    end = time.time()

    mem_after = process.memory_info().rss
    peak_mem = max(mem_before, mem_after) / 1024 ** 3  # in GB
    delta_mem = (mem_after - mem_before) / 1024 ** 3

    safe_log(
        f"Peak memory usage: {peak_mem:.2f} GB (Δ {delta_mem:+.2f} GB), time: {end - start:.2f}s"
    )
    return result

def reindex_adata_to_reference(adata: ad.AnnData, reference_var_df: pd.DataFrame) -> ad.AnnData:
    """
    Reindex an AnnData object to align with a reference .var DataFrame.

    This ensures that:
    - All genes in the reference appear in the returned object (with zeros for missing ones).
    - The `.var` matches the reference metadata exactly (by index and content).

    Parameters
    ----------
    adata : ad.AnnData
        The input AnnData object (already filtered/processed).
    reference_var_df : pd.DataFrame
        A full .var DataFrame including all expected genes (with metadata).

    Returns
    -------
    ad.AnnData
        A new AnnData object with padded expression matrix and full reference .var.
    """
    #safe_log("Reindexing AnnData to reference var_names (optimized)")

    reference_var_names = reference_var_df.index
    adata_genes = adata.var_names

    missing = reference_var_names.difference(adata_genes)
    present = reference_var_names.intersection(adata_genes)

    # Extract only the present genes in the correct order
    X_present = adata[:, present].X

    # Create zero matrix for missing genes
    if len(missing) > 0:
        X_missing = csr_matrix((adata.n_obs, len(missing)))
        X_final = hstack([X_present, X_missing], format="csr")
        del X_present, X_missing
    else:
        X_final = X_present

    # Reorder columns to match reference
    combined_gene_order = list(present) + list(missing)
    reorder_idx = pd.Index(combined_gene_order).get_indexer(reference_var_names)
    X_final = X_final[:, reorder_idx]
    del reorder_idx, combined_gene_order

    # Sanity check
    assert X_final.shape[1] == len(reference_var_names), "Mismatch in final gene count"

    # Ensure the .var DataFrame has a proper index name
    var_out = reference_var_df.copy(deep=True)
    if var_out.index.name is None:
        var_out.index.name = "gene_id"

    return ad.AnnData(
        X=X_final,
        obs=adata.obs.copy(),
        var=var_out
    )

def batch_filter_write_zarr(
    h5ad_files: List[str],
    expressions: List[str],
    group_cols: List[str],
    count_df: pd.DataFrame,
    batch_size: int,
    output_dir: str,
    log_memory_flag: bool = False,
    use_raw_counts: bool = True
) -> List[str]:
    """
        Filter and batch-process .h5ad files into Zarr chunks with harmonized gene space.

        This function:
        - Computes the union of all gene metadata (.var) across .h5ad files.
        - Applies cell filters (if any) to each AnnData.
        - Groups data into batches, reindexes to include all genes (zero-padding missing ones).
        - Writes Zarr chunks to disk for later merging or MTX export.

        Parameters
        ----------
        h5ad_files : List[str]
            Paths of h5ad files to process.
        expressions : List[str]
            List of Python expressions to filter cells (e.g., 'adata.obs["n_genes"] > 500').
        group_cols : List[str]
            List of column names in .obs to use for grouping (e.g., ['donor_id']).
        count_df : pd.DataFrame
            DataFrame with counts of cells per group for downsampling.
        batch_size : int
            Number of filtered AnnData objects to group per Zarr chunk.
        output_dir : str
            Output directory to write Zarr chunks (one per batch).
        log_memory_flag : bool
            If True, logs memory usage after each batch is processed.
        use_raw_counts : bool
            If True and `.raw` exists, uses raw counts matrix for filtering and writing.

        Returns
        -------
        List[str]
            List of paths to Zarr chunks that were written.

        Notes
        -----
        - Reindexing ensures that every output chunk has the same .var structure and same gene ordering.
        - All Zarr chunks will share the same reference_var_df (gene metadata), even if input files are heterogeneous.
        - Batches that contain no matching cells after filtering are skipped.
        - Memory usage is logged if `log_memory_flag` is enabled.

        Raises
        ------
        SystemExit
            If any filter expression fails to evaluate or yields incorrect shape.
        AssertionError
            If reindexing produces mismatched var_names compared to the reference.
        """

    reference_var_df = get_union_var_df(h5ad_files)
    logger.info(f"Total unique genes across all files: {reference_var_df.shape[0]}")

    os.makedirs(output_dir, exist_ok=True)

    chunk_paths = []
    batch = []
    batch_index = 0

    logger.info("Gathering data from files in batches")

    def finalize_batch(batch: list, index: int, log_memory_flag:bool) -> str:
        combined = ad.concat(batch, join='outer', merge='same')
        if log_memory_flag:
            combined = measure_peak_memory(reindex_adata_to_reference, combined, reference_var_df)
        else:
            combined = reindex_adata_to_reference(combined, reference_var_df)

        # sanitize after reindexing
        combined.obs = sanitize_dataframe(combined.obs)
        combined.var = sanitize_dataframe(combined.var)
        # safety check
        assert list(combined.var_names) == list(reference_var_df.index), "Final batch var_names mismatch"

        if log_memory_flag:
            log_memory(f"after preparing batch [{index + 1}] with [{combined.n_obs}] cells")

        path = os.path.join(output_dir, f"chunk_{index}.zarr")
        combined.write_zarr(store=path, convert_strings_to_categoricals=False)
        return path

    pbar = tqdm(h5ad_files, total=len(h5ad_files), desc="Normalizing and filtering inputs")

    for fname in pbar:
        adata = ad.read_h5ad(fname)

        if use_raw_counts and adata.raw is not None:
            adata = extract_raw_counts(adata, cast_to_int=True)

        adata_subset = filters.evaluate_filters(adata, expressions)

        if adata_subset is not None and count_df is not None and len(group_cols)>0:
            adata_subset = downsample_anndata(
                adata=adata_subset,
                group_cols=group_cols,
                count_df=count_df,
                filename=os.path.basename(fname),  # Match what's in your count_df
                seed=42
            )
        if adata_subset is not None:
            batch.append(adata_subset)
        else:
            logger.info(f"{fname}: no cells matched filters")

        if len(batch) >= batch_size:
            chunk_path = finalize_batch(batch, batch_index, log_memory_flag)
            chunk_paths.append(chunk_path)
            batch.clear()
            batch_index += 1

    if batch:
        chunk_path = finalize_batch(batch, batch_index, log_memory_flag)
        chunk_paths.append(chunk_path)

    return chunk_paths

def write_10x_mtx_from_zarr_chunks(
    zarr_chunk_paths: List[str],
    output_dir: str,
    var_reference_column: Optional[str] = None,
    ten_x_format: bool = True,
) -> None:
    """
    Write a 10x-style MTX (Matrix Market) directory by streaming from Zarr chunks.

    Parameters
    ----------
    zarr_chunk_paths : list of str
        Paths to Zarr chunk directories that have aligned `.var_names`.
    output_dir : str
        Directory to write matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz.
    var_reference_column : str, optional
        If provided, use this column from `.var` for gene_id and gene_name
        instead of the index. Default is None.
    """
    os.makedirs(output_dir, exist_ok=True)

    matrix_path = os.path.join(output_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(output_dir, "barcodes.tsv.gz")
    features_path = os.path.join(output_dir, "features.tsv.gz")

    # -------------------------------
    # Pass 1: Gather metadata
    # -------------------------------
    n_cells = 0
    nnz = 0
    all_barcodes = []
    reference_var_names = None
    reference_var_df = None

    for path in tqdm(zarr_chunk_paths, desc="Pass 1: Gathering cell/gene/counts metadata"):
        adata = ad.read_zarr(path)
        nnz += adata.X.nnz
        n_cells += adata.n_obs
        all_barcodes.extend(adata.obs_names)

        if reference_var_names is None:
            reference_var_names = list(adata.var_names)
            reference_var_df = adata.var
        else:
            if list(adata.var_names) != reference_var_names:
                raise ValueError(f"Chunk {path} has different var_names than the reference")

    n_genes = len(reference_var_names)

    # Write barcodes.tsv.gz
    pd.Series(all_barcodes).to_csv(
        barcodes_path, index=False, header=False, sep="\t", compression="gzip"
    )

    # Write features.tsv.gz
    if var_reference_column:
        if reference_var_df is None or var_reference_column not in reference_var_df.columns:
            raise ValueError(f"Column {var_reference_column!r} not found in reference .var DataFrame")
        # gene_name comes from the requested column…
        gene_names = reference_var_df[var_reference_column].tolist()
        # …but gene_id is always the DataFrame’s index
        gene_ids = reference_var_df.index.astype(str).tolist()
    else:
        # no reference DataFrame → both name and id fall back to whatever you passed in
        gene_names = reference_var_names
        gene_ids = reference_var_names

    features_df = pd.DataFrame({
        "gene_id": gene_ids,
        "gene_name": gene_names,
        "feature_type": "Gene Expression"
    })

    features_df.to_csv(
        features_path, index=False, header=False, sep="\t", compression="gzip"
    )

    # -------------------------------
    # Pass 2: Stream matrix.mtx.gz
    # -------------------------------
    if ten_x_format:
        write_10x_matrix_market_file_faster(
            matrix_path=matrix_path,
            zarr_chunk_paths=zarr_chunk_paths,
            n_cells=n_cells,
            n_genes=n_genes,
            nnz=nnz,
            flush_every=100_000
        )
    else:
        write_standard_matrix_market_file(
            matrix_path=matrix_path,
            zarr_chunk_paths=zarr_chunk_paths,
            n_cells=n_cells,
            n_genes=n_genes,
            nnz=nnz,
            flush_every=100_000
        )

def sparse_equal(A, B):
    A = A.tocsr() if issparse(A) else csr_matrix(A)
    B = B.tocsr() if issparse(B) else csr_matrix(B)
    return (A != B).nnz == 0

def read_mtx_native_format(mtx_dir: str) -> ad.AnnData:
    """
    Read a Matrix Market MTX file (native orientation: cells × genes) and construct an AnnData object.

    Assumes the following files exist in the directory:
    - matrix.mtx or matrix.mtx.gz
    - barcodes.tsv or barcodes.tsv.gz
    - features.tsv or features.tsv.gz

    The features.tsv file must contain exactly 3 columns: gene_id, gene_name, feature_type.

    Parameters
    ----------
    mtx_dir : str
        Path to directory containing the matrix.mtx(.gz), barcodes.tsv(.gz), and features.tsv(.gz) files.

    Returns
    -------
    AnnData
        An AnnData object with `.X` as sparse matrix, `.obs_names` from barcodes, `.var` from features.
    """
    def read_tsv(filename: str) -> pd.DataFrame:
        path = os.path.join(mtx_dir, filename)
        compression = "gzip" if path.endswith(".gz") else None
        return pd.read_csv(path, sep="\t", header=None, compression=compression)

    # Read matrix
    matrix_file = next(f for f in os.listdir(mtx_dir) if f.startswith("matrix.mtx"))
    X = io.mmread(os.path.join(mtx_dir, matrix_file)).tocsr()

    # Read barcodes
    barcodes_file = next(f for f in os.listdir(mtx_dir) if f.startswith("barcodes.tsv"))
    barcodes = read_tsv(barcodes_file)[0].astype(str)

    # Read features
    features_file = next(f for f in os.listdir(mtx_dir) if f.startswith("features.tsv"))
    features_df = read_tsv(features_file)
    if features_df.shape[1] != 3:
        raise ValueError("features.tsv must have exactly 3 columns: gene_id, gene_name, feature_type")
    features_df.columns = ["gene_id", "gene_name", "feature_type"]
    features_df.index = features_df["gene_id"]

    # Construct AnnData
    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes.tolist()
    adata.var = features_df
    adata.var_names = features_df.index.tolist()

    return adata


def collect_unique_obs_var_columns(h5ad_files) -> tuple:
    """
    Read each H5AD file in backed='r' mode, and collect unique obs and var column names.

    Parameters
    ----------
    h5ad_files : list of str
        Paths to .h5ad files.
    expressions : list of str, optional

    Returns
    -------
    tuple of (set, set)
        A tuple containing the set of all obs column names and the set of all var column names.
    """
    unique_obs_cols = set()
    unique_var_cols = set()

    for fp in tqdm(h5ad_files, desc="Loading H5AD files"):
        adata = ad.read_h5ad(fp, backed="r")
        unique_obs_cols.update(adata.obs.columns)
        unique_var_cols.update(adata.var.columns)
        # close the backing file to free resources
        try:
            adata.file.close()
        except Exception:
            pass

    return unique_obs_cols, unique_var_cols


def view_column_values(h5ad_files, column_spec, expressions):
    """
    Print a summary of values for a single obs or var column across all h5ad_files.

    column_spec must be "obs.<colname>" or "var.<colname>".
    expressions is a list of Python filter expressions (only applied to obs).
    """
    try:
        axis, colname = column_spec.split(".", 1)
    except ValueError:
        raise ValueError("`--view_values` must be prefixed with 'obs.<col>' or 'var.<col>'")

    if axis not in ("obs", "var"):
        raise ValueError("Prefix must be 'obs' or 'var' (e.g. obs.biobank or var.gene_id)")

    # Only show filters if we're doing obs
    if axis == "obs":
        if expressions:
            print("Applying filter expressions to obs:")
            for expr in expressions:
                print(f"  - {expr}")
        else:
            print("No filter expressions provided; using full obs DataFrame.")
    else:
        # var mode: ignore expressions
        if expressions:
            print("Warning: --expr filters are ignored when viewing var columns.")
        print("Gathering var values without filtering.")

    all_arrays = []
    for fp in tqdm(h5ad_files, desc=f"Gathering {column_spec}"):
        adata_raw = ad.read_h5ad(fp, backed="r")

        # pick the right DataFrame and optionally filter
        try:
            if axis == "obs":
                df = adata_raw.obs
                df_subset = filters.filter_obs_dataframe(df, expressions)
            else:
                df_subset = adata_raw.var
        except ValueError as e:
            # Friendly exit on bad filter (including missing columns)
            print(f"\nError in filter: {e}")
            sys.exit(1)

        if colname not in df_subset.columns:
            print(f"[!] Column {colname!r} not found in {axis} of {os.path.basename(fp)}")
            try: adata_raw.file.close()
            except: pass
            continue

        all_arrays.append(df_subset[colname].values)
        try: adata_raw.file.close()
        except: pass

    if not all_arrays:
        print("No values found for", column_spec)
        return

    flat = np.concatenate(all_arrays)
    series = pd.Series(flat)

    print(f"\nSummary for {column_spec} ({len(flat)} values):")
    if pd.api.types.is_numeric_dtype(series):
        print(f"  min:    {series.min()}")
        print(f"  max:    {series.max()}")
        print(f"  mean:   {series.mean()}")
        print(f"  median: {series.median()}")
    else:
        counts = series.value_counts()
        print("  value counts:")
        for val, cnt in counts.items():
            print(f"    {val!r}: {cnt}")

def main() -> None:
    """
    Entry point for command-line tool. Parses arguments, filters AnnData, and writes output.
    """
    parser = argparse.ArgumentParser(description="Filter and merge AnnData files based on user-defined expressions.")
    parser.add_argument("--dir", required=True, help="Directory containing .h5ad files.")
    parser.add_argument("--expr", action='append', default=[],
                        help="Optional Python expressions to filter cells (e.g., 'adata.obs[\"n_genes\"] > 500').")

    parser.add_argument("--group-cols", action="append", required=False,
                        help="Repeatable: column(s) in .obs to use for grouping (e.g. --group-cols donor_id)."
                             "This is used to downsample cells to a maximum number per group across all files.")
    parser.add_argument("--max-cells", type=int, required=False,
                        help="Maximum number of cells to retain per group across all files.")
    parser.add_argument("--downsample_report", help="File containing the counts of cells per group for "
                                                    "downsampling.", required=False)

    parser.add_argument("--output", required=False, help="Output file (.h5ad) or directory (MTX).")
    parser.add_argument("--output-format", choices=["h5ad", "10x_mtx", "mtx"], default="h5ad",
                        help="Output format: h5ad or mtx (10x-style) or mtx (h5ad native shape format). Default is h5ad.")

    parser.add_argument("--mtx-var-name", default=None,
                        help="If the output format is MTX, this is the name of the column in .var to use for gene_id "
                             "and gene_name. Default is None, which uses the index.")

    parser.add_argument("--use-raw", action="store_true",
                        help="Use raw counts instead of normalized data.")
    parser.add_argument("--batch-size", type=int, default=5,
                        help="Number of files to combine per chunk when writing temp Zarr files. Default is 5. " +
                             "More chunks increase memory usage and speed.")
    parser.add_argument("--log-level", default="INFO", help="Logging level (e.g., DEBUG, INFO, WARNING)")
    parser.add_argument("--log-memory", action="store_true",
                        help="If set, logs memory usage after major processing steps.")

    parser.add_argument("--view_columns", action="store_true", default=False, help="If set, will print the unique column names "
                                                                                   "of the .obs and .var DataFrames from all anndata inputs, then exit.  "
                                                                                   "All other arguments are ignored.")
    parser.add_argument("--view_values", required=False, help="If set, will print a summary of the obs or var column values across all files, then exit."
                                                              " All other arguments are ignored.  Prefix the column name with 'obs.' or 'var,' to specify which DataFrame to view.")


    args = parser.parse_args()

    try:
        args = parser.parse_args(sys.argv[1:])
    except SystemExit as e:
        print("Argparse failed with error:")
        parser.print_help()
        raise

    if not any("ipykernel" in m for m in sys.modules):
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s] %(message)s",
                level=logging.INFO,
                datefmt='%Y-%m-%d %H:%M:%S'
            )

    # validate args
    if args.group_cols is not None and args.max_cells is None:
        raise ValueError("If group-cols is specified, max-cells must also be specified.")

    # Get the input files
    h5ad_files = sorted(glob.glob(os.path.join(args.dir, "*.h5ad")))

    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad files found in directory: {args.dir}")

    #if in view mode, all other arguments are ignored and the program quits after printing the columns.
    if args.view_columns:
        obs_cols, var_cols = collect_unique_obs_var_columns(h5ad_files)

        print("\nUnique obs columns:")
        for col in sorted(obs_cols):
            print(col)

        print("\nUnique var columns:")
        for col in sorted(var_cols):
            print(col)

        sys.exit(0)

    #if in view column value mode, all other arguments are ignored and the program quits after printing the values.
    if args.view_values:
        view_column_values(h5ad_files, args.view_values, args.expr)
        return 0

    #Not in view mode, check the output file is set
    if args.output is None:
        parser.error("An output file or directory must be specified unless in view mode.")

    # if downsampling is requested, compute it, otherwise it's None.
    rate_df = None
    if args.group_cols is not None and args.max_cells is not None:
        logger.info(f"Downsampling to {args.max_cells} cell barcodes per group requested.  "
                    f"Calculating exact downsampling requirements")
        rate_df = compute_downsampling_counts(
            h5ad_files=h5ad_files,
            group_cols=args.group_cols,
            max_cells=args.max_cells,
            expressions=args.expr
        )
        if args.downsample_report:
            logger.info(f"Writing downsampling report to {args.downsample_report}")
            rate_df.to_csv(args.downsample_report, sep="\t", index=False)

    logger.info(f"Filtering AnnData files in directory: {args.dir}")

    result=None
    with tempfile.TemporaryDirectory() as tmpdir:
        batch_zarr_dir = os.path.join(tmpdir, "zarr_chunks")
        logger.info(f"Using temporary Zarr directory: {batch_zarr_dir}")
        chunk_paths = batch_filter_write_zarr(
            h5ad_files,
            args.expr,
            args.group_cols,
            rate_df,
            batch_size=args.batch_size,
            output_dir=batch_zarr_dir,
            log_memory_flag=args.log_memory,
            use_raw_counts=args.use_raw
        )

        if args.log_memory:
            log_memory("after writing Zarr chunks")
        # merge if MTX
        if args.output_format == "10x_mtx":
            write_10x_mtx_from_zarr_chunks(chunk_paths, args.output, args.mtx_var_name, ten_x_format=True)
            logger.info(f"Using 10x format: Saved matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz to {args.output}")
        elif args.output_format == "mtx":
            write_10x_mtx_from_zarr_chunks(chunk_paths, args.output, args.mtx_var_name, ten_x_format=False)
            logger.info(f"Using mtx format: Saved matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz to {args.output}")

        # don't merge if output format is mtx
        elif args.output_format == "h5ad":
            #slightly smarter concatonation where you only read one chunk at a time.
            result = ad.read_zarr(chunk_paths[0])
            #t1 = as_dataframe(result)
            for path in chunk_paths[1:]:
                next_chunk = ad.read_zarr(path)
                result = ad.concat([result, next_chunk], join='outer', merge='same')

            #every time I concat I'm sanitizing the obs.
            result.obs = sanitize_dataframe(result.obs)
            result.obs = sanitize_dataframe(result.obs)
            logger.info(f"Final result shape: {result.shape}, X type: {type(result.X)}")
            if args.log_memory:
                log_memory("after final concat")
            result.write_h5ad(args.output)
        else:
            raise ValueError(f"Unknown output format: {args.output_format}")

    logger.info("Finished processing all files.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
