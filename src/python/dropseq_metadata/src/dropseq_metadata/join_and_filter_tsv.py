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
Load a primary tab-separated file, and zero or more secondary tab-separated files, each with join columns
specified.  Each secondary file is left outer joined to the primary file on the join column.

The values in the join column of a join file must be unique.

If there is a collision between the column names in the primary and secondary files, the column in the primary
file is used.  If there is a collision between the column names in the secondary files, the column in the earlier
file is used.

In addition, there may be zero or more (column name, column value) pairs that are set on each line unconditionally.

Each resulting row is passed through zero more filters, which can apply a min or a max threshold to a column,
have a list of values to match against, or have a list of values to exclude, which may be specified on the command
line or in a file.
"""

import argparse
import sys
import pandas as pd

DELETEME_COLUMN_SUFFIX = '_deleteme'

def try_convert_string_to_number(s):
    int_val = None
    float_val = None
    try:
        int_val = int(s)
    except ValueError:
        pass
    try:
        float_val = float(s)
    except ValueError:
        pass
    if int_val is None and float_val is None:
        return s
    if int_val is None and float_val is not None:
        return float_val
    if float_val is None and int_val is not None:
        return int_val
    if int_val == float_val:
        return int_val
    return float_val

def load_values_file(file):
    return pd.read_csv(file, sep='\t', header=None).iloc[0]

def add_subparser(subparsers):
    parser = subparsers.add_parser("join_and_filter_tsv", description=__doc__)
    parser.add_argument("--input", "-i", type=argparse.FileType('r'),
                        help="Primary tab-separated file to join.  Default: %(default)s", default=sys.stdin)
    parser.add_argument("--output", "-o", type=argparse.FileType('w'),
                        help="Output file.  Default: %(default)s", default=sys.stdout)
    parser.add_argument("--join", "-j", nargs=3, action='append', default=[],
                        metavar=('SECONDARY_FILE', 'INPUT_COL', 'JOIN_COL'),
                        help="Secondary tab-separated file to join, and the columns in the primary input "
                             "and join file to join on. May be specified multiple times.")
    parser.add_argument("--set", "-s", nargs=2, action='append', default=[], metavar=('COLUMN', 'VALUE'),
                        help="Set a column to a constant value.  May be specified multiple times.")
    parser.add_argument("--min", nargs=2, action='append', default=[], metavar=('COLUMN', 'VALUE'),
                        help="Filter out rows where COLUMN is less than VALUE.  May be specified multiple times.")
    parser.add_argument("--max", nargs=2, action='append', default=[], metavar=('COLUMN', 'VALUE'),
                        help="Filter out rows where COLUMN is greater than VALUE.  May be specified multiple times.")
    parser.add_argument("--include-file", nargs=2, action='append', default=[], metavar=('COLUMN', 'FILE'),
                        help="Filter out rows where COLUMN is not in FILE.  May be specified multiple times.")
    parser.add_argument("--exclude-file", nargs=2, action='append', default=[], metavar=('COLUMN', 'FILE'),
                        help="Filter out rows where COLUMN is in FILE.  May be specified multiple times.")
    parser.add_argument("--include", nargs='+', action='append', default=[], metavar=('COLUMN', 'VALUE'),
                        help="Filter out rows where COLUMN is not one of the given VALUEs.  May be specified multiple times.")
    parser.add_argument("--exclude", nargs='+', action='append', default=[], metavar=('COLUMN', 'VALUE'),
                        help="Filter out rows where COLUMN is one of the given VALUEs.  May be specified multiple times.")

def main(options):
    # load the primary file
    primary = pd.read_csv(options.input, sep='\t')
    options.input.close()
    # load each secondary file, and join it to the primary file, dropping secondary columns that are already in the primary
    for join_file, input_col, join_col in options.join:
        join_col_in_left = join_col in primary.columns
        secondary = pd.read_csv(join_file, sep='\t')
        primary = primary.merge(secondary, how='left', left_on=input_col, right_on=join_col,
                                suffixes=(None, DELETEME_COLUMN_SUFFIX))
        if not join_col_in_left:
            # drop the join column from the merged data frame
            primary.drop(join_col, axis=1, inplace=True)
        # drop the secondary columns that are already in the primary
        for col in primary.columns:
            if col.endswith(DELETEME_COLUMN_SUFFIX):
                primary.drop(col, axis=1, inplace=True)
    # set columns to constant values
    for column, value in options.set:
        primary[column] = try_convert_string_to_number(value)
    # filter out rows based on column values
    for column, value in options.min:
        primary = primary[primary[column] >= try_convert_string_to_number(value)]
    for column, value in options.max:
        primary = primary[primary[column] <= try_convert_string_to_number(value)]
    for column, file in options.include_file:
        include_values = load_values_file(file)
        primary = primary[primary[column].isin(include_values)]
    for column, file in options.exclude_file:
        exclude_values = load_values_file(file)
        primary = primary[~primary[column].isin(exclude_values)]
    for includes in options.include:
        column = includes[0]
        values = includes[1:]
        values = [try_convert_string_to_number(value) for value in values]
        primary = primary[primary[column].isin(values)]
    for excludes in options.exclude:
        column = excludes[0]
        values = excludes[1:]
        values = [try_convert_string_to_number(value) for value in values]
        primary = primary[~primary[column].isin(values)]
    # write the output
    primary.to_csv(options.output, sep='\t', index=False)
    options.output.close()
    return 0


