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
"""
Load a dense tabular DGE, optionally convert gene names to Ensembl gene IDs, and write to h5ad file.
"""

import argparse
import anndata
import scanpy as sc
import logging
from scipy import sparse
import pandas as pd
import sys


def read_dropseq_dge(in_file: str, dtype: str = "uint") -> sc.AnnData:
    # copied from CellAnnotationService.ipynb
    logging.info(f"Reading in dropseq DGE file {in_file}")
    dge = sc.read_text(in_file, delimiter="\t", dtype=dtype)
    dge = dge.transpose()
    dge.X = sparse.csr_matrix(dge.X)
    return dge


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input", "-i", required=True,
                        help="Tabular DGE file, optionally gzipped")
    parser.add_argument("--output", "-o", required=True, help="Output h5ad file")
    parser.add_argument("--reduced-gtf", "-r", help="Tabular file with columns gene_id and gene_name "
                                                    "(and any other columns).  If provided, gene names will be "
                                                    "converted to gene IDs")
    parser.add_argument("--dtype", "-d", default="uint32", help="Data type for the matrix.  Default: %(default)s)")
    options = parser.parse_args(args)

    dge = read_dropseq_dge(options.input, dtype=options.dtype)
    if options.reduced_gtf:
        logging.info(f"Converting gene names to gene IDs using {options.reduced_gtf}")
        gtf = pd.read_csv(options.reduced_gtf, sep="\t")
        # filter gtf to only include annotationType="gene" in order to avoid multiple entries for the same gene name
        gtf = gtf[gtf['annotationType'] == "gene"]
        # find gene names for which there is not a gene ID
        missing_gene_ids = set(dge.var_names) - set(gtf['gene_name'])
        if missing_gene_ids:
            raise Exception(f"Missing gene IDs for {len(missing_gene_ids)} genes: " + ", ".join(missing_gene_ids))
        dge.var["gene_symbol"] = dge.var_names
        dge.var_names = dge.var_names.map(gtf.set_index('gene_name')['gene_id'])

    logging.info(f"Writing {options.output}")
    adDge = anndata.AnnData(X=dge.X, obs=dge.obs, var=dge.var, dtype=options.dtype)
    adDge.write_h5ad(options.output, compression="gzip")
    return 0


if __name__ == "__main__":
    sys.exit(main())
