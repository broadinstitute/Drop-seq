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
Simple python tools for converting between hdf5 and tabular digital gene expression file formats.
"""
import argparse
import logging
import sys

try:
    from . import hdf5_10X_to_text
    from . import optimus_h5ad_to_dropseq
except ImportError:
    import hdf5_10X_to_text
    import optimus_h5ad_to_dropseq



def main(args=None):
    parser = argparse.ArgumentParser(prog="dropseq_hdf5", description=__doc__)
    parser.add_argument("--log-level", "-l", default="INFO", choices=dctLogLevel.keys(),
                        help="Set the logging level.  (default: %(default)s)")
    subparsers = parser.add_subparsers(
        title="sub-commands",
        description="valid commands",
        dest="tool")
    hdf5_10X_to_text.add_subparser(subparsers)
    optimus_h5ad_to_dropseq.add_subparser(subparsers)

    if args is None:
        args = sys.argv[1:]
    if len(args) == 0:
        parser.print_help()
        return 1
    else:
        options = parser.parse_args(args)
        logger.setLevel(dctLogLevel[options.log_level])
        if options.tool == "hdf5_10X_to_text":
            return hdf5_10X_to_text.main(options)
        elif options.tool == "optimus_h5ad_to_dropseq":
            return optimus_h5ad_to_dropseq.main(options)
        else:
            # should be unpossible because parse_args will complain
            raise ValueError(f"Unrecognized tool: {options.tool}")


if __name__ == "__main__":
    sys.exit(main())
