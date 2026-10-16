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

Any of the options can also be read from a file with --args-file, which has one argument per line, taken exactly as
written (no quoting, no trimming; a blank line is an empty argument).  The arguments of the file replace the
--args-file option where it appears, so the result is the same as typing them there, and the order of repeated options
is kept.  --args-file may be given more than once, but a file cannot name another --args-file.
"""

import argparse
import sys

import pandas as pd
from pandas.errors import MergeError
from dropseq.aggregation import logger, add_log_argument

DELETEME_COLUMN_SUFFIX = '_deleteme'

BOOLEAN_STRINGS = {'true': True,
                     'false': False}

def try_convert_string(s):
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
        if s.lower() in BOOLEAN_STRINGS:
            return BOOLEAN_STRINGS[s.lower()]
        return s
    if int_val is None and float_val is not None:
        return float_val
    if float_val is None and int_val is not None:
        return int_val
    if int_val == float_val:
        return int_val
    return float_val

def load_values_file(file):
    return pd.read_csv(file, sep='\t', header=None).iloc[:,0]

ARGS_FILE_OPTION = "--args-file"


def read_args_file(parser, file_name):
    """
    :return: list of the arguments in the file, one per line, verbatim
    """
    try:
        with open(file_name, newline='') as f:
            text = f.read()
    except OSError as e:
        parser.error(f"Cannot read {ARGS_FILE_OPTION} {file_name}: {e}")
    lines = text.split('\n')
    if lines[-1] == '':
        # the newline that ends the last line (or an empty file)
        lines.pop()
    for line in lines:
        if line == ARGS_FILE_OPTION or line.startswith(ARGS_FILE_OPTION + '='):
            parser.error(f"{ARGS_FILE_OPTION} {file_name} cannot contain {ARGS_FILE_OPTION}")
    return lines


def expand_args_files(parser, args):
    """
    :return: args with each --args-file FILE (or --args-file=FILE) replaced by the arguments in the file
    """
    expanded = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == ARGS_FILE_OPTION:
            if i + 1 >= len(args):
                parser.error(f"argument {ARGS_FILE_OPTION}: expected one argument")
            expanded += read_args_file(parser, args[i + 1])
            i += 2
        elif arg.startswith(ARGS_FILE_OPTION + '='):
            expanded += read_args_file(parser, arg[len(ARGS_FILE_OPTION) + 1:])
            i += 1
        else:
            expanded.append(arg)
            i += 1
    return expanded


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
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
    parser.add_argument("--set-first", action='store_true', default=False,
                        help="Put the columns specified with --set first, in the order given, rather than last.")
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
    parser.add_argument("--drop", action='append', default=[],
                        help="Column to drop from the output.  May be specified multiple times.")
    parser.add_argument("--rename", nargs=2, action='append', default=[], metavar=('OLD', 'NEW'),
                        help="Rename column.  May be specified multiple times.")
    parser.add_argument("--no-header", action='store_true', default=False,
                        help="Do not write a header row.")
    # Expanded before parsing, so this only documents the option (and catches an abbreviation, see below).
    parser.add_argument(ARGS_FILE_OPTION, action='append', default=[], metavar='FILE',
                        help="File of arguments, one per line, taken verbatim and used in place of this option.  "
                             "May be specified multiple times.  The file cannot contain --args-file.")
    options = parser.parse_args(expand_args_files(parser, sys.argv[1:] if args is None else args))
    if options.args_file:
        # only reachable through an abbreviation such as --args-f, which would otherwise be silently ignored
        parser.error(f"write {ARGS_FILE_OPTION} in full")
    del options.args_file
    return options

def main(args=None):
    run(parse_args(args))

def run(options):
    # load the primary file
    primary = pd.read_csv(options.input, sep='\t', comment='#')
    options.input.close()
    # load each secondary file, and join it to the primary file, dropping secondary columns that are already in the primary
    for join_file, input_col, join_col in options.join:
        join_col_in_left = join_col in primary.columns
        secondary = pd.read_csv(join_file, sep='\t', comment='#')
        try:
            # many_to_one: require that the join column in the secondary file is unique
            primary = primary.merge(secondary, how='left', left_on=input_col, right_on=join_col, validate="many_to_one",
                                    suffixes=(None, DELETEME_COLUMN_SUFFIX))
        except MergeError    as e:
            logger.error(f"Error joining {join_file} on {input_col} and {join_col}: {e}")
            return 1
        if not join_col_in_left:
            # drop the join column from the merged data frame
            primary.drop(join_col, axis=1, inplace=True)
        # drop the secondary columns that are already in the primary
        for col in primary.columns:
            if col.endswith(DELETEME_COLUMN_SUFFIX):
                primary.drop(col, axis=1, inplace=True)
    # set columns to constant values
    for column, value in options.set:
        primary[column] = try_convert_string(value)
    # filter out rows based on column values
    for column, value in options.min:
        primary = primary[primary[column] >= try_convert_string(value)]
    for column, value in options.max:
        primary = primary[primary[column] <= try_convert_string(value)]
    for column, file in options.include_file:
        include_values = load_values_file(file)
        primary = primary[primary[column].isin(include_values)]
    for column, file in options.exclude_file:
        exclude_values = load_values_file(file)
        primary = primary[~primary[column].isin(exclude_values)]
    for includes in options.include:
        column = includes[0]
        values = includes[1:]
        values = [try_convert_string(value) for value in values]
        primary = primary[primary[column].isin(values)]
    for excludes in options.exclude:
        column = excludes[0]
        values = excludes[1:]
        values = [try_convert_string(value) for value in values]
        primary = primary[~primary[column].isin(values)]
    # drop columns
    for column in options.drop:
        primary = primary.drop(column, axis=1)
    # rename columns
    for column, value in options.rename:
        primary.rename(columns={column: value}, inplace=True)
    # move the columns set to constant values to the front
    if options.set_first:
        set_columns = list(dict.fromkeys(column for column, _ in options.set if column in primary.columns))
        primary = primary[set_columns + [c for c in primary.columns if c not in set_columns]]
    # write the output
    primary.to_csv(options.output, sep='\t', index=False, header=not options.no_header)
    options.output.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
