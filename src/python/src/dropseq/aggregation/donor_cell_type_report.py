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
Report the number and fraction of nuclei, and the median UMIs per nucleus, of each cell type of each donor.

This is a port of DropSeq.aggregation::donorCellTypeReport, with snake_case column names (donor, cell_type,
num_nuclei, fraction_nuclei, median_umis_per_nucleus).  Cells with no donor are ignored.  Cells with a donor but no cell
type, for example because they have no MapMyCells result, count toward the number of cells of the donor, which is the
denominator of fraction_nuclei, but have no row.  (The R function counted such cells as one extra nucleus of an NA cell
type.)  Donors and cell types are in order of first appearance in the input.
"""

import argparse
import sys

import pandas as pd

from dropseq.aggregation import logger, add_log_argument

OUTPUT_COLUMNS = ["donor", "cell_type", "num_nuclei", "fraction_nuclei", "median_umis_per_nucleus"]


def donor_cell_type_report(cell_metadata, cell_type_column, donor_column="donor",
                           umi_column="num_retained_transcripts"):
    """
    :param cell_metadata: DataFrame of cell metadata
    :return: DataFrame with OUTPUT_COLUMNS, one row per donor and cell type
    """
    cells = cell_metadata.dropna(subset=[donor_column])
    rows = []
    for donor, donor_cells in cells.groupby(donor_column, sort=False):
        typed = donor_cells.dropna(subset=[cell_type_column])
        for cell_type, type_cells in typed.groupby(cell_type_column, sort=False):
            rows.append((donor, cell_type, len(type_cells), round(len(type_cells) / len(donor_cells), 3),
                         type_cells[umi_column].median()))
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    parser.add_argument("--input", "-i", required=True, help="Tab-separated cell metadata file.")
    parser.add_argument("--output", "-o", required=True, help="Tab-separated report.")
    parser.add_argument("--cell-type-column", required=True, help="Cell metadata column with the cell type.")
    parser.add_argument("--donor-column", default="donor", help="Cell metadata column with the donor.")
    parser.add_argument("--umi-column", default="num_retained_transcripts",
                        help="Cell metadata column with the number of UMIs.")
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    # donors and cell types are labels, even when they look like numbers
    columns = [options.donor_column, options.cell_type_column]
    cell_metadata = pd.read_csv(options.input, sep="\t", dtype={c: str for c in columns})
    missing = [c for c in [*columns, options.umi_column] if c not in cell_metadata.columns]
    if missing:
        logger.error(f"No {missing} column(s) in {options.input}")
        return 1
    report = donor_cell_type_report(cell_metadata, options.cell_type_column, options.donor_column,
                                    options.umi_column)
    report.to_csv(options.output, sep="\t", index=False, float_format="%.15g")
    return 0


if __name__ == "__main__":
    sys.exit(main())
