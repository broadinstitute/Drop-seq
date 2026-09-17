#!/usr/bin/env python3
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
Read a MapMyCells csv and count the number of cells of each supercluster_name.
Write a tab-separated file with the counts, and fraction of cells to the output file.
"""

import argparse
import sys

import pandas as pd


def count_mmc_celltypes(input_csv) -> pd.DataFrame:
    """
    Count the number of cells of each supercluster_name in a MapMyCells csv file.

    Parameters
    ----------
    input_csv
        Path or file-like object for the MapMyCells csv file.

    Returns
    -------
    pd.DataFrame
        A data frame indexed by supercluster_name with 'count' and 'fraction' columns,
        sorted by count in descending order.
    """
    df = pd.read_csv(input_csv, comment='#', dtype=str)
    counts = df['supercluster_name'].value_counts()
    counts.name = 'count'
    result = counts.to_frame()
    result['fraction'] = (result['count'] / result['count'].sum()).round(4)
    return result


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', '-i', default=sys.stdin, type=argparse.FileType('r'),
                        help='MapMyCells csv file.  Default: stdin.')
    parser.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'),
                        help='Output tab-separated file.  Default: stdout.')
    options = parser.parse_args(args)

    result = count_mmc_celltypes(options.input)
    result.to_csv(options.output, sep='\t', index_label='supercluster_name')
    options.output.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

