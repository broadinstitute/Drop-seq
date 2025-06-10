#!/usr/bin/env python3
# MIT License
#
# Copyright 2022 Broad Institute
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
Merge parquet files together to either another parquet file or a tab-delimited file.
"""
import argparse
import sys
from types import MethodType

import pandas as pd

import misc_utils.pandas_utils as pd_utils
from misc_utils.argparse_utils import argparse_error


def main(args=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        allow_abbrev=False,
    )
    parser.error = MethodType(argparse_error, parser)

    parser.add_argument(
        '--input',
        '-i',
        action='append',
        required=True,
        help='Input parquet files.',
    )
    parser.add_argument(
        '--output',
        '-o',
        default=sys.stdout.buffer,
        type=argparse.FileType('wb'),
        help='Output file. If not specified, output to stdout.',
    )
    parser.add_argument(
        '--output-format',
        '-f',
        choices=['parquet', 'tsv'],
        default=None,
        help='Output format. If not specified, output tab separated values.',
    )

    options = parser.parse_args(args)
    dfs = [pd.read_parquet(fIn) for fIn in options.input]
    df = pd.concat(dfs)

    output_format = options.output_format
    if not output_format:
        if options.output == sys.stdout.buffer:
            output_format = 'tsv'
        elif options.output.name.endswith('.parquet'):
            output_format = 'parquet'
        else:
            output_format = 'tsv'

    if output_format == 'parquet':
        df.to_parquet(options.output)
    else:
        pd_utils.to_tsv(df, options.output, index=False)


if __name__ == "__main__":
    sys.exit(main())
