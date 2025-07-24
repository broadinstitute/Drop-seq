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
# Annotate eQTLs with information from the GTF, dbSNP, and eventually other sources.
"""

import argparse
import sys
from types import MethodType
import pandas as pd
import dropseq.metadata.read_gtf as read_gtf
import dropseq.hdf5.io_utils as io_utils

from dropseq.util.argparse_utils import argparse_error
import dropseq.util.pandas_utils as pd_utils


class Columns:
    gene_id = "gene_id"
    phenotype_id = "phenotype_id"

def join_gtf(qtls, gtfPath, missing_value):
    """
    Annotate eQTLs with gene information from a GTF file.
    """
    # Read the GTF file, filtering for 'gene' features only
    gtf = read_gtf.read_gtf(gtfPath, lambda df: df[df[read_gtf.GtfRequiredColNames.FEATURE] == 'gene'])
    # select relevant columns
    gtf = gtf[gtf.FEATURE == 'gene']
    gtf = gtf[['gene_name', read_gtf.GtfRequiredColNames.CHROMOSOME, read_gtf.GtfRequiredColNames.START,
               read_gtf.GtfRequiredColNames.END, read_gtf.GtfRequiredColNames.STRAND, 'gene_id']]
    # rename columns to match what we want to go into the qtls DataFrame
    gtf.rename(columns={
        read_gtf.GtfRequiredColNames.CHROMOSOME: 'gene_chr',
        read_gtf.GtfRequiredColNames.START: 'gene_start',
        read_gtf.GtfRequiredColNames.END: 'gene_end',
        read_gtf.GtfRequiredColNames.STRAND: 'strand',
    }, inplace=True)
    # change these columns to strings so that missing values can be replaced with a string
    gtf['gene_start'] = gtf['gene_start'].astype(str)
    gtf['gene_end'] = gtf['gene_end'].astype(str)
    gtf['strand'] = gtf['strand'].astype(str)
    qtls = qtls.merge(gtf, on=Columns.gene_id, how='left')
    qtls.fillna(missing_value, inplace=True)
    return qtls

def join_dbsnp(qtls, dbsnpPath, missing_value):
    colnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    chromIndex = colnames.index('CHROM')
    posIndex = colnames.index('POS')
    idIndex = colnames.index('ID')
    refIndex = colnames.index('REF')
    altIndex = colnames.index('ALT')
    numColsOfInterest = len(colnames)
    # slurp dbSNP file, parse important columns, and put into a dictionary with variant_id constructed from dbSNP columns.
    dctDbSnp = {}
    with open(dbsnpPath, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith('#'):
                continue
            fields = line.split('\t', maxsplit=numColsOfInterest)
            if len(fields) < numColsOfInterest:
                raise ValueError(f"Error parsing dbSNP file {dbsnpPath}.  Not enough column in line '{line}'.")
            alternateAlleles = fields[altIndex].split(',')
            for alternateAllele in alternateAlleles:
                variant_id = ":".join([fields[chromIndex], fields[posIndex], fields[refIndex], alternateAllele])
                dctDbSnp[variant_id] = fields[idIndex]
    # add rsids to qtls DataFrame
    rsids = [dctDbSnp.get(variant_id, missing_value) for variant_id in qtls['variant_id']]
    qtls['rsid'] = rsids
    return qtls



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
        required=True,
        help='Input cis_qtl.txt.gz file produced by tensorQTL.',
    )
    parser.add_argument(
        '--output',
        '-o',
        help='Output file. If not specified, output to stdout.',
    )

    parser.add_argument(
        '--gtf',
        '-g',
        help='GTF file to annotate eQTLs with gene information.',
    )

    parser.add_argument(
        '--dbsnp',
        '-d',
        help='dbSNP file to annotate eQTLs with SNP information.',
    )
    parser.add_argument('--missing-value', '-m', default='.',
                        help='String to use for missing values, typically rsids, gene_name, gene_start, gene_end. (default: %(default)s)')
    options = parser.parse_args(args)
    if options.output is None:
        fOut = sys.stdout
    else:
        fOut = io_utils.open_maybe_gz(options.output, 'wt')
    qtls = pd.read_csv(options.input, sep='\t')
    qtls.rename(columns={
        Columns.phenotype_id: Columns.gene_id
    }, inplace=True)
    if options.gtf:
        qtls = join_gtf(qtls, options.gtf, options.missing_value)
    if options.dbsnp:
        qtls = join_dbsnp(qtls, options.dbsnp, options.missing_value)

    if options.output is None:
        fOut = sys.stdout
    else:
        fOut = io_utils.open_maybe_gz(options.output, 'wt')
    for k, v in vars(options).items():
        if k != 'output' and v is not None:
            fOut.write(f"# {k}={v}\n")
    qtls.to_csv(fOut, sep='\t', index=False)
    if options.output is not None:
        fOut.close()
    return 0



if __name__ == "__main__":
    sys.exit(main())
