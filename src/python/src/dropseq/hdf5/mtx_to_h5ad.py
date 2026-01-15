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
Convert matrix.mtx(.gz) to h5ad file.
"""

import argparse
import logging
import sys

import scanpy as sc


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input-directory", "-i", required=True,
                        help="Directory containing matrix.mtx[.gz], barcodes.tsv[.gz], and features.tsv[.gz]")
    parser.add_argument("--output", "-o", required=True, help="Output h5ad file")
    parser.add_argument("--prefix", "-p",
                        help="Optional prefix for matrix, barcodes, and features files.  Default: none")
    options = parser.parse_args(args)

    adDge = sc.read_10x_mtx(options.input_directory, prefix=options.prefix)
    # seems like downstream tools want csr format rather than csc.  It would be nice if there were a standard.
    adDge.X = adDge.X.tocsr()

    logging.info(f"Writing {options.output}")
    adDge.write_h5ad(options.output, compression="gzip")
    return 0


if __name__ == "__main__":
    sys.exit(main())
