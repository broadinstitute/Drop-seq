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
Convert a MapMyCells csv to a tab-separated file without the '#' header lines, optionally prefixing the column
names (except cell_id, which is the join key) with a string, typically '<MMC model>_'.

This is a port of DropSeq.cellclassification::convertPrefixMMCCsvClp.  The output can be joined to a cell metadata
file with join_and_filter_tsv.

--cell-type-column-output writes the name of the cell type column, <prefix><level>_name, after checking that it is in
the output.  The level is --level, or the top level of the taxonomy hierarchy of the csv.
"""

import argparse
import json
import re
import sys

import pandas as pd

from dropseq.aggregation import logger, add_log_argument

CELL_ID_COLUMN = "cell_id"
TAXONOMY_HIERARCHY_PATTERN = re.compile(r"^#\s*taxonomy hierarchy\s*=\s*(\[.*\])\s*$")
# R writes numbers with up to 15 significant digits, so 1.0000 is written as 1.
FLOAT_FORMAT = "%.15g"


def read_taxonomy_hierarchy(mmc_csv):
    """
    Parse the '# taxonomy hierarchy = ["neighborhood", "class", ...]' header line of a MapMyCells csv.

    :param mmc_csv: path of the MapMyCells csv.
    :return: list of taxonomy levels, from the top level down.
    """
    with open(mmc_csv) as f:
        for line in f:
            if not line.startswith("#"):
                break
            match = TAXONOMY_HIERARCHY_PATTERN.match(line.strip())
            if match:
                return json.loads(match.group(1))
    raise ValueError(f"No '# taxonomy hierarchy = [...]' header line found in {mmc_csv}")


def convert_mmc_csv(mmc_csv, column_prefix=None, columns_to_skip=(CELL_ID_COLUMN,)):
    """
    :param mmc_csv: path or file-like object of the MapMyCells csv.
    :param column_prefix: prefix for the column names, or None for no prefix.
    :param columns_to_skip: columns that are not prefixed.
    :return: DataFrame
    """
    df = pd.read_csv(mmc_csv, comment="#")
    if column_prefix:
        df.columns = [c if c in columns_to_skip else column_prefix + c for c in df.columns]
    return df


def cell_type_column(mmc_csv, df, column_prefix=None, level=None):
    """
    :param mmc_csv: path of the MapMyCells csv, for its taxonomy hierarchy.
    :param df: the converted DataFrame.
    :param level: taxonomy level.  Default: the top level of the taxonomy hierarchy.
    :return: name of the cell type column of df, <prefix><level>_name
    """
    if level is None:
        level = read_taxonomy_hierarchy(mmc_csv)[0]
    column = f"{column_prefix or ''}{level}_name"
    if column not in df.columns:
        raise ValueError(f"No {column} column in {mmc_csv}")
    return column


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    parser.add_argument("--input", "-i", required=True, help="MapMyCells csv file.")
    parser.add_argument("--output", "-o", type=argparse.FileType("w"), default=sys.stdout,
                        help="Output tab-separated file.  Default: stdout.")
    parser.add_argument("--column-prefix", "-p", default=None,
                        help="Prefix for every column name except cell_id, e.g. HMBA_Human_WB_v0.5_.")
    parser.add_argument("--level", default=None,
                        help="Taxonomy level of the cell type.  Default: the top level of the taxonomy hierarchy.")
    parser.add_argument("--cell-type-column-output", default=None,
                        help="Write the name of the cell type column, <prefix><level>_name, to this file.")
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    df = convert_mmc_csv(options.input, options.column_prefix)
    if options.cell_type_column_output:
        try:
            column = cell_type_column(options.input, df, options.column_prefix, options.level)
        except ValueError as e:
            logger.error(str(e))
            return 1
        with open(options.cell_type_column_output, "w") as f:
            f.write(column + "\n")
    df.to_csv(options.output, sep="\t", index=False, float_format=FLOAT_FORMAT)
    options.output.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
