#!/usr/bin/env python3
# MIT License
# 
# Copyright 2025 Broad Institute
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
Swap the gene IDs var["gene_symbol"] and names (var_names) in an h5ad file.

This can be useful when running MapMyCells with query markers that use gene names.
"""
import argparse
import sys
import anndata as ad

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input", "-i", required=True, help="Input h5ad file")
    parser.add_argument("--output", "-o", required=True, help="Output h5ad file")
    options = parser.parse_args(args)
    adata = ad.read_h5ad(options.input)
    gene_symbol_colname = "gene_symbol"
    if gene_symbol_colname not in adata.var.columns:
        raise ValueError(f"Input h5ad file does not contain '{gene_symbol_colname}' in var columns.")
    adata.var['gene_id'] = adata.var_names.to_series()
    adata.var_names = adata.var["gene_symbol"]
    adata.write_h5ad(options.output, compression="gzip")
    return 0


if __name__ == "__main__":
    sys.exit(main())
