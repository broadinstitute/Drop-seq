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
Return 0 if the given h5ad, features.tsx.gz, tabular DGE or MapMyCells query markers json file has Ensembl IDs.
"""

import argparse
import gzip
import json
import re
import sys

import anndata as ad

def open_maybe_gz(file_path):
    """
    Open a file, handling .gz files automatically.
    :param file_path: Path to the file
    :return: File object
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def isEnsemblId(gene):
    """
    Check if the given gene ID is an Ensembl ID.
    Ensembl IDs typically start with 'ENSG' followed by 11 digits.
    :param gene: Gene ID to check
    :return: True if the gene ID is an Ensembl ID, False otherwise
    """
    return bool(re.fullmatch(r"ENSG\d{11}", gene))

def getGeneFromFeaturesFile(features_path):
    """
    Read the features file and return the first gene ID (first column).
    :param features_path: Path to the features.tsv(.gz) file
    :return: string from the first column of the first row
    """
    with open_maybe_gz(features_path) as f:
        line = f.readline().rstrip("\n")
    parts = line.split("\t")
    if len(parts) != 3:
        raise Exception(f"features file {features_path} should have 3 columns, but has {len(parts)}: {line}")
    return parts[0]

def getGeneFromH5ad(h5ad_path):
    """
    Read the h5ad file and return the first gene ID from the var_names.
    :param h5ad_path: Path to the h5ad file
    :return: string from the first row of var_names
    """
    adata = ad.read_h5ad(h5ad_path, backed='r')
    return adata.var_names[0]

def getGeneFromDgeFile(dge_path):
    with open_maybe_gz(dge_path) as f:
        seenHeader = False
        for line in f:
            if not line.startswith("#"):
                parts = line.rstrip("\n").split("\t")
                if seenHeader:
                    return parts[0]
                elif parts[0] != "GENE":
                    raise Exception(f"Expected first line of DGE file {dge_path} to be a header line starting with 'GENE', but got: {line}")
                else:
                    seenHeader = True
    raise Exception(f"Got to end of DGE file {dge_path} without seeing a gene.")

def getGeneFromQueryMarkers(query_markers_path):
    """
    Read the MapMyCells query markers json file and return the first gene ID.
    :param query_markers_path: Path to the MapMyCells query markers json file
    :return: string from the first gene in the 'genes' list
    """
    with open(query_markers_path) as f:
        data = json.load(f)
    for markerList in data.values():
        if len(markerList) > 0:
            return markerList[0]
    raise Exception(f"Query markers file {query_markers_path} does not contain any non-empty list.")

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--h5ad", "-i", help="Input h5ad file")
    group.add_argument("--features", "-f", help="Input features.tsv.gz file")
    group.add_argument("--dge", "-d", help="Input tabular DGE file")
    group.add_argument("--query-markers", "-q", help="Input MapMyCells query markers json file")
    options = parser.parse_args(args)
    if options.h5ad:
        gene = getGeneFromH5ad(options.h5ad)
    elif options.features:
        gene = getGeneFromFeaturesFile(options.features)
    elif options.dge:
        gene = getGeneFromDgeFile(options.dge)
    elif options.query_markers:
        gene = getGeneFromQueryMarkers(options.query_markers)
    else:
        raise ValueError("No input file specified")
    if isEnsemblId(gene):
        print(f"Gene ID {gene} is an Ensembl ID.", file=sys.stderr)
        return 0
    else:
        print(f"Gene ID {gene} is NOT an Ensembl ID.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
