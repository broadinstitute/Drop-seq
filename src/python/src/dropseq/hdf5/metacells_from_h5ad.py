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
import argparse
import os
import logging
import sys

import anndata as ad
import scanpy as sc
import pandas as pd

from scipy.sparse import issparse
from tqdm import tqdm

import dropseq.hdf5.filters as filters


logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Create metacells for an input h5ad.")
    parser.add_argument("--input", required=True, type=str, help="File with paths to h5ads to include.")
    parser.add_argument("--output", required=True, type=str, help="Output metacell h5ad.")
    parser.add_argument("--groups", required=True, type=str, nargs="+", help="Columns to create metacells for.")
    parser.add_argument("--filtering", action='append', default=[], help="Filtering to apply to h5ad.")
    parser.add_argument("--gene-coding", required=False, default="gene_symbol", type=str, help="How to represent genes in the output, must exist in h5ad var.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")

    return parser.parse_args()

def validate_inputs(args):
    """
    Validate command-line arguments.
    
    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments. 
    
    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    """
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file {args.input} not found.")
    else:
        if args.verbose: 
            logger.info(f"Input file {args.input} found.")

def read_h5ad_list(input: str) -> list[str]:
    """
    Read a list of file paths from a text file.
    
    Parameters
    ----------
    input : str
        Path to the text file containing the list of file paths.

    Returns
    -------
    list[str]
        List of file paths.
    """
    with open(input, 'r') as file:
        file_list = [line.strip() for line in file]

    return(file_list)

def read_single_h5ad(h5ad_path: str, verbose: bool) -> ad.AnnData:
    """
    Read an AnnData object from a file.
    
    Parameters
    ----------
    h5ad_path : str
        Path to the h5ad file.
    verbose : bool
        If True, log the reading process.
    
    Returns
    -------
    ad.AnnData
        The AnnData object read from the file."""

    if verbose: 
        logger.info(f"Reading AnnData object from {h5ad_path}")
    adata = sc.read(h5ad_path)
    return adata

def filter_adata(adata: ad.AnnData, filtering: list[str], verbose: bool) -> ad.AnnData:
    """
    Filter the AnnData object based on the specified filtering criteria (if provided).
    If no filtering criteria is provided, the original AnnData object is returned.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to filter.
    filtering : list[str]
        List of filtering criteria to apply to the AnnData object.
    verbose : bool
        If True, log the filtering process.
    Returns
    -------
    ad.AnnData
        The filtered AnnData object.
    """
    if filtering is None:
        if verbose: 
            logger.info("No filtering criteria provided.")
        return adata
    else: 
        if verbose:
            logger.info(f"Applying filtering criteria: {filtering}")
        filtered_adata = filters.evaluate_any_filters(adata, filtering)
        return filtered_adata
    
def validate_groups(adata: ad.AnnData, groups: list[str], verbose: bool) -> None:
    """
    Validate that the specified groups exist in the AnnData object.
    
    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to validate.
    groups : list[str]
        List of groups to check for in the AnnData object.
    verbose : bool
        If True, log the validation process.
    
    Raises
    ------
    ValueError
        If any of the specified groups are not found in the AnnData object.
    """  
    for group in groups:
        if group not in adata.obs.columns:
            raise ValueError(f"Group {group} not found in adata.obs.")
        elif verbose:
                logger.info(f"Group {group} found in adata.obs.")

def reindex_adata_var(adata: ad.AnnData, index: str, verbose: bool) -> ad.AnnData:
    """
    Reindex the AnnData object based on the specified index.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to reindex.
    index : str
        The column name to use as the new index.

    Returns
    -------
    ad.AnnData
        The reindexed AnnData object.
    """
    
    if adata.var.index.name == index:
        if verbose:
            logger.info(f"Index {index} already set as index.")
    elif index not in adata.var.columns:
        raise ValueError(f"Index {index} not found in adata.var.columns.")
    else:
        if verbose:
            logger.info(f"Reindexing AnnData object using {index} as index.")
        adata.var_names = adata.var[index]    
    return adata

def reformat_adata_groups(adata: ad.AnnData, groups: list[str]) -> ad.AnnData:
    """
    Replace spaces with hyphens in the grouping columns of the AnnData object.

    This improves the output format for the metacells, such that an output column
    would look like "MD6366_STR-D1-MSN-GABA" instead of "MD6366_STR D1 MSN GABA".

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to reformat.
    groups : list[str]
        List of groups to reformat in the AnnData object.
    verbose : bool
        If True, log the reformatting process.

    Returns
    -------
    ad.AnnData
        The reformatted AnnData object.
    """

    adata.obs[groups]=adata.obs[groups].apply(lambda col: col.str.replace(' ', '-', regex=False) 
                                              if col.dtype == 'category' else col)
    return adata

def create_metacells(adata: ad.AnnData, groups: list[str], verbose: bool) -> ad.AnnData:
    """
    Create metacells by aggregating the AnnData object based on the specified groups.
    
    This uses `scanpy.get.aggregate` which will return the metacells as an AnnData object.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to aggregate.
    groups : list[str]
        List of groups to aggregate by.
    verbose : bool
        If True, log the aggregation process.
    
    Returns
    -------
    ad.AnnData
        The aggregated AnnData object (metacells).
    """
    if verbose:
        logger.info(f"Creating metacells using groups: {groups}")
    metacells = sc.get.aggregate(adata, by=groups, func=["sum"])
    metacells.X = metacells.layers["sum"]
    if verbose:
        logger.info(f"Metacells shape: {metacells.shape}")
    return metacells

def concatenate_metacells(metacells_list: list[ad.AnnData], verbose: bool) -> ad.AnnData:
    """
    Concatenate a list of AnnData objects into a single AnnData object using an outer join
    and replace NA values with zeros.

    Parameters:
    -----------
    metacells_list : list[ad.AnnData]
        List of AnnData objects to concatenate.
    verbose : bool
        If True, log the concatenation process.

    Returns:
    --------
    ad.AnnData
        The concatenated AnnData object.
    """

    if len(metacells_list) == 1:
        if verbose: 
            logger.info("Only one metacells object found, returning it.")
        return metacells_list[0]
    else:
        if verbose:
            logger.info(f"Concatenating {len(metacells_list)} metacells objects.")
        concatenated = ad.concat(metacells_list, axis=0, join="outer")
        concatenated.X = pd.DataFrame(concatenated.X).fillna(0).values
        return concatenated


def convert_adata_to_dataframe(adata: ad.AnnData, transpose: bool=False) -> pd.DataFrame:
    """
    Convert an AnnData object to a pandas DataFrame.
    
    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to convert.
    transpose : bool, optional
        If True, transpose the DataFrame (default is False).
    
    Returns
    -------
    pd.DataFrame
        The converted DataFrame.
    """

    if issparse(adata.X):
        df = pd.DataFrame.sparse.from_spmatrix(
            adata.X,
            index=adata.obs_names, 
            columns=adata.var_names
        )
    else:
        df = pd.DataFrame(
            adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        )

    if transpose:
        return df.T
    else:
        return df

def collapse_metacells(metacells: ad.AnnData) -> ad.AnnData:
    """
    Collapse the metacells by summing the values of duplicate rows.
    
    Parameters
    ----------
    metacells : ad.AnnData
        The AnnData object containing the metacells.
    
    Returns
    -------
    ad.AnnData
        The AnnData object with collapsed metacells.
    """

    metacell_df = convert_adata_to_dataframe(metacells)
    metacell_df_summed = metacell_df.groupby(metacell_df.index).sum()
    metacell_df_summed = metacell_df_summed.astype(int)
    metacell_adata = ad.AnnData(X=metacell_df_summed, obs=pd.DataFrame(index=metacell_df_summed.index), var=metacells.var.copy())
    metacell_adata = metacell_adata[metacell_adata.obs_names.sort_values()]
    metacell_adata.var = metacell_adata.var.sort_index()
    return metacell_adata

def write_h5ad(metacells: ad.AnnData, output: str) -> None:
    """Write the metacells to an h5ad file.
    
    Parameters
    ----------
    metacells : ad.AnnData
        The AnnData object containing the metacells.
    output : str
        The output file path. Must end with .h5ad.

    Returns
    -------
    None
    """
    logger.info(f"Writing metacells to {output}")
    metacells.write_h5ad(output)

def write_txt(metacells: ad.AnnData, output: str) -> None:
    """Write the metacells to a txt file.
    
    Parameters
    ----------
    metacells : ad.AnnData  
        The AnnData object containing the metacells.
    output : str
        The output file path. Must end with .txt or .txt.gz.

    """
    df = convert_adata_to_dataframe(metacells, transpose=True)

    if output.endswith(".txt.gz"):
        logger.info(f"Writing metacells to {output}")
        df.to_csv(output, sep='\t', compression='gzip', index=True)
    elif output.endswith(".txt"):
        logger.info(f"Writing metacells to {output}")
        df.to_csv(output, sep='\t', index=True)
    else:
        logger.error("Output file name does not end with .txt or .txt.gz")

def write_metacells(metacells: ad.AnnData, output: str) -> None:
    """
    Write the metacells to a file.
    
    The output file format is determined by the file extension. 
    .h5ad, .txt, or .txt.gz are supported.

    Parameters
    ----------
    metacells : ad.AnnData
        The AnnData object containing the metacells.
    output : str
        The output file path. Must end with .h5ad, .txt, or .txt.gz.
    
    Raises
    ------
    ValueError
        If the output file does not end with .h5ad, .txt, or .txt.gz.

    """
    if output.endswith(".h5ad"):
        write_h5ad(metacells, output)
    elif output.endswith(".txt") or output.endswith(".txt.gz"):
        write_txt(metacells, output)
    else:
        raise ValueError("Output file must end with .h5ad, .txt, or .txt.gz")

def generate_metacells(h5ad_list: list[str], 
                       filtering: list[str],
                       groups: list[str],
                       gene_coding: str,
                       verbose: bool) -> ad.AnnData:
    """
    Main workflow function, generates metacells from a list of h5ad files.
    
    Parameters
    ----------
    h5ad_list : list[str]
        List of paths to h5ad files.
    filtering : list[str]
        List of filtering criteria to apply to the h5ad files.
    groups : list[str]  
        List of groups to create metacells for.
    gene_coding : str 
        How to represent genes in the output, must exist in h5ad var (e.g. gene_symbol).
    verbose : bool
        If True, log the process.
    
    Returns
    -------
    ad.AnnData
        The final metacells AnnData object."""

    metacells_list = []

    for h5ad_path in tqdm(h5ad_list, desc="Generating metacells from individual h5ads"):
        logger.info(f"Processing {h5ad_path}")
        adata = read_single_h5ad(h5ad_path, verbose)
        adata = filter_adata(adata, filtering, verbose)
        validate_groups(adata, groups, verbose)
        adata=reindex_adata_var(adata, gene_coding, verbose)
        adata=reformat_adata_groups(adata, groups)
        metacells = create_metacells(adata, groups, verbose)
        metacells_list.append(metacells)

    if verbose:
        logger.info(f"Number of metacells to aggregate: {len(metacells_list)}")
        logger.info("Concatenating metacells...")
    
    metacells = concatenate_metacells(metacells_list, verbose)
    
    if verbose:
        logger.info("Metacells concatenated.")
        logger.info("Collapsing duplicate rows in metacells...")
    
    metacells = collapse_metacells(metacells)

    if verbose:
        logger.info("Duplicate rows summed.")
    
    logger.info(f"Final metacells shape: {metacells.shape}")

    return metacells

class TqdmLoggingHandler(logging.StreamHandler):
    def __init__(self, stream=None):
        super().__init__(stream)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logging_with_tqdm():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if not any(isinstance(h, TqdmLoggingHandler) for h in logger.handlers):
        handler = TqdmLoggingHandler()
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

def main():
    setup_logging_with_tqdm()
    args = parse_arguments()
    validate_inputs(args)
    h5ad_list = read_h5ad_list(args.input)
    metacells = generate_metacells(h5ad_list, args.filtering, args.groups, args.gene_coding, args.verbose)
    write_metacells(metacells, args.output)
    logger.info("Metacells generated successfully.")

if __name__ == "__main__":
    main()

