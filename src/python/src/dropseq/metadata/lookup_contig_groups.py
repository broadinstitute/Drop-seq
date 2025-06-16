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
Retrieve one or more contig groups from a contig_groups.yaml by looking up the groups in the labels.
"""

import argparse
import sys

import yaml

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--contig-groups',
        '-c',
        help='A yaml file containing labels for all the contigs in the reference to be created.  '
             'Use the program CreateContigGroups to make a skeleton yaml, then edit it according to the '
             'instructions in the file.',
        required=True,
    )
    parser.add_argument(
        '--group',
        '-g',
        help='One or more contig groups',
        action='extend',
        nargs='+',
        default=[],
    )
    parser.add_argument(
        '--out',
        '-o',
        default='-',
        type=argparse.FileType('w', encoding='UTF-8'),
        help='where to write the output',
    )
    args = parser.parse_args()
    return(run(args))

def run(args):
    try:
        with open(args.contig_groups) as fIn:
            contig_groups = yaml.load(fIn, Loader=yaml.SafeLoader)
        for group in args.group:
            for contig, labels in contig_groups.items():
                if not isinstance(labels, list):
                    labels = [labels]
                if group in labels:
                    args.out.write(contig + '\n')
    finally:
        args.out.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
