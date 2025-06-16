#!/usr/bin/env python3
# MIT License
#
# Copyright 2020 Broad Institute
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
Convert DGE matrix in 10X hdf5 format produced by remove-background into dense Drop-seq DGE.
"""
import argparse
import gzip

import numpy

try:
    from . import io_utils
    from . import downstream
    from . import cli
except ImportError:
    import io_utils
    import downstream
    import cli



def readDgeHeader(dge_path):
    retval = []
    if dge_path is not None:
        with io_utils.open_maybe_gz(dge_path) as fIn:
            for strLine in fIn:
                if not strLine.startswith("#"):
                    break
                retval.append(strLine.rstrip("\n"))
    return retval


def readCbrbCommandLineFromLog(cbrb_log):
    sawPrecedingLine = False
    with open(cbrb_log) as fIn:
        for strLine in fIn:
            strLine = strLine.rstrip("\n")
            if sawPrecedingLine:
                return strLine
            elif strLine.endswith("Command:"):
                sawPrecedingLine = True
    raise Exception("Could not find command line in " + cbrb_log)

def add_subparser(subparsers):
    parser = subparsers.add_parser("hdf5_10X_to_text", description=__doc__)
    parser.add_argument("--input", "-i", required=True, help="Input hdf5 file.")
    parser.add_argument("--output", "-o", required=True,
                        help="Output tab-separated text file, with rows=genes and columns=cells.  "
                             "If .gz extension, will be gzipped.")
    parser.add_argument("--output-sizes", type=argparse.FileType(mode="w"),
                        help="If specified, output a tab-separated text file with columns CELL_BARCODE and NUM_TRANSCRIPTS")
    parser.add_argument("--progress-interval", "-p", default=1000, type=int,
                        help="Report progress after this many genes.  Set to zero to suppress progress messages.  (default: %(default)s)")
    parser.add_argument("--analyzed-barcodes-only", default=False, action='store_true',
                        help="output a limited set of barcodes: only those analyzed by the algorithm.  (default: output all barcodes)")
    parser.add_argument("--limit", "-l", type=int, help="Output no more than this number of genes (for debugging)")
    parser.add_argument("--header", help="If set, read DGE header lines from this file and write to output.")
    parser.add_argument("--cbrb-log",
                                         help="If set, read the CBRB log for the CBRB "
                                              "command line, and write to #COMMAND record in DGE header.")


def main(options):
    h5 = downstream.anndata_from_h5(options.input, analyzed_barcodes_only=options.analyzed_barcodes_only)
    genes = [g for g in h5.var_names]
    cells = [c for c in h5.obs_names]
    mat = h5.X  # This has cell rows and gene columns
    mat = mat.astype(numpy.int32)
    # Convert to column-oriented because we want to access by column.
    mat = mat.tocsc()

    if options.output.endswith(".gz"):
        fOut = gzip.open(options.output, mode="wt", encoding="ascii")
    else:
        fOut = open(options.output, "w")

    header_lines = readDgeHeader(options.header)
    for header_line in header_lines:
        print(header_line, file=fOut)

    cbrb_command_line = None
    if options.cbrb_log is not None:
        cbrb_command_line = readCbrbCommandLineFromLog(options.cbrb_log)
    if cbrb_command_line is not None:
        print("#COMMAND\tCL:" + cbrb_command_line, file=fOut)

    if options.output_sizes is not None:
        cell_sizes = [0] * len(cells)
    cli.logger.info("%d genes to write" % len(genes))
    print("\t".join(["GENE"] + cells), file=fOut)
    num_genes = len(genes)
    if options.limit is not None and options.limit < num_genes:
        num_genes = options.limit
    for i in range(num_genes):
        expressionForGene = mat[:, i].toarray().flatten().tolist()
        if options.output_sizes is not None:
            for j in range(len(cell_sizes)):
                cell_sizes[j] += expressionForGene[j]
        print("\t".join([genes[i]] + [str(v) for v in expressionForGene]), file=fOut)
        if options.progress_interval > 0 and (i + 1) % options.progress_interval == 0:
            cli.logger.info("Wrote %d genes" % (i + 1))
    fOut.close()

    if options.output_sizes is not None:
        print("\t".join(["cell_barcode", "num_transcripts"]), file=options.output_sizes)
        for j in range(len(cell_sizes)):
            print("\t".join([cells[j], str(cell_sizes[j])]), file=options.output_sizes)
        options.output_sizes.close()
    cli.logger.info("Done")
    return 0
