# MIT License
#
# Copyright 2026 Broad Institute
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
"""
Build the h5ad file of one library from its DGE and its cell metadata.

X holds the counts of the DGE (cells x genes, sparse) for the cells in the cell metadata, in the order of the cell
metadata.  obs is the cell metadata, indexed by <PREFIX>_<cell barcode>, the cell names that MakeTripletDge gives the
cells of the aggregated matrix, so that the h5ad files of several libraries can be combined.  var has the column gene_symbol, with
the gene symbols of the DGE, and is indexed by gene ID (index name gene_ids) if there is a reduced GTF, otherwise by
gene symbol.  This is the layout of the h5ad files of the CAP freeze.  obs and var are sanitized with
aggregate_filter_adata.sanitize_dataframe, so that missing values are encoded in a way that h5ad can store.
"""

import argparse
import gzip
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import dropseq.hdf5.h5ad_defs as h5ad_defs
from dropseq.aggregation import logger, add_log_argument
from dropseq.hdf5.aggregate_filter_adata import sanitize_dataframe

# labels, even when they look like numbers
GENE_SYMBOL_COLUMN = "gene_symbol"
STRING_COLUMNS = ["PREFIX", "cell_barcode", "donor"]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_dge(path, barcodes):
    """
    Read the columns of a dense DGE for some cells, one gene at a time so that the dense matrix is never in memory.

    :param path: Drop-seq DGE, optionally gzipped, with leading '#' header lines, then a header line starting with
    GENE and the cell barcodes, then a line per gene
    :param barcodes: cell barcodes to keep, in the order of the result
    :return: genes (list of gene names), CSR matrix of cells (the barcodes) x genes
    """
    genes = []
    rows, cols, values = [], [], []
    with open_text(path) as dge:
        line = dge.readline()
        while line.startswith("#"):
            line = dge.readline()
        header = line.rstrip("\n").split("\t")
        position = {barcode: i for i, barcode in enumerate(header[1:])}
        missing = [b for b in barcodes if b not in position]
        if missing:
            raise ValueError(f"{len(missing)} of the {len(barcodes)} cells of the cell metadata are not in {path}, "
                             f"for example {missing[:3]}")
        keep = np.array([position[b] for b in barcodes]) + 1  # +1 for the gene column
        for line in dge:
            fields = line.rstrip("\n").split("\t")
            counts = np.array(fields, dtype=object)[keep].astype(float)
            nonzero = np.flatnonzero(counts)
            rows.extend([len(genes)] * len(nonzero))
            cols.extend(nonzero)
            values.extend(counts[nonzero])
            genes.append(fields[0])
    # genes x cells, then transposed to the cells x genes of h5ad
    matrix = sp.coo_matrix((values, (rows, cols)), shape=(len(genes), len(barcodes))).T.tocsr()
    if np.all(np.equal(np.mod(matrix.data, 1), 0)):
        matrix = matrix.astype(np.uint32)
    return genes, matrix


def make_var(gene_names, reduced_gtf=None):
    """
    :param gene_names: gene symbols of the DGE
    :param reduced_gtf: optional reduced GTF, with columns gene_name, gene_id and annotationType, to find gene IDs
    :return: var DataFrame with the column gene_symbol, indexed by gene ID (named gene_ids) if there is a reduced GTF,
    otherwise by gene symbol
    """
    var = pd.DataFrame({GENE_SYMBOL_COLUMN: list(gene_names)}, index=pd.Index(gene_names))
    if reduced_gtf:
        gtf = pd.read_csv(reduced_gtf, sep="\t", usecols=["gene_name", "gene_id", "annotationType"])
        gtf = gtf[gtf["annotationType"] == "gene"].drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
        missing = [g for g in gene_names if g not in gtf.index]
        if missing:
            raise ValueError(f"{len(missing)} genes of the DGE are not in {reduced_gtf}, for example {missing[:3]}")
        var.index = pd.Index(gtf.loc[list(gene_names)].to_numpy(), name=h5ad_defs.GENE_IDS_COLUMN_NAME)
    if not var.index.is_unique:
        raise ValueError(f"Duplicate {'gene IDs' if reduced_gtf else 'gene names'}, for example "
                         f"{var.index[var.index.duplicated()][:3].tolist()}")
    return var


def build_library_h5ad(dge, cell_metadata, reduced_gtf=None, prefix_column="PREFIX", barcode_column="cell_barcode"):
    """
    :param cell_metadata: DataFrame with one row per cell
    :return: AnnData with obs from cell_metadata, indexed by <prefix>_<barcode>, and with the genes of the DGE
    """
    obs = cell_metadata.copy()
    if not obs[barcode_column].is_unique:
        raise ValueError(f"Duplicate {barcode_column} in the cell metadata")
    obs.index = pd.Index((obs[prefix_column] + "_" + obs[barcode_column]).to_numpy())
    genes, matrix = read_dge(dge, obs[barcode_column].tolist())
    var = make_var(genes, reduced_gtf)
    return ad.AnnData(X=matrix, obs=sanitize_dataframe(obs), var=sanitize_dataframe(var))


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    parser.add_argument("--dge", required=True, help="Tabular DGE of the library, optionally gzipped.")
    parser.add_argument("--cell-metadata", required=True,
                        help="Tab-separated cell metadata of the library, one row per cell.")
    parser.add_argument("--reduced-gtf",
                        help="Reduced GTF, with columns gene_name, gene_id and annotationType.  If given, var is "
                             "indexed by gene ID (the ENSG IDs) and has the gene symbols in the gene_symbol column.  "
                             "Otherwise var is indexed by the gene symbols of the DGE.")
    parser.add_argument("--output", "-o", required=True, help="Output h5ad file.")
    parser.add_argument("--prefix-column", default="PREFIX", help="Cell metadata column with the library name.")
    parser.add_argument("--barcode-column", default="cell_barcode", help="Cell metadata column with the cell barcode.")
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    dtype = {c: str for c in [*STRING_COLUMNS, options.prefix_column, options.barcode_column]}
    cell_metadata = pd.read_csv(options.cell_metadata, sep="\t", dtype=dtype)
    missing = [c for c in [options.prefix_column, options.barcode_column] if c not in cell_metadata.columns]
    if missing:
        logger.error(f"No {missing} column(s) in {options.cell_metadata}")
        return 1
    adata = build_library_h5ad(options.dge, cell_metadata, options.reduced_gtf, options.prefix_column, options.barcode_column)
    logger.info(f"Writing {adata.n_obs} cells x {adata.n_vars} genes to {options.output}")
    adata.write_h5ad(options.output, compression="gzip")
    return 0


if __name__ == "__main__":
    sys.exit(main())
