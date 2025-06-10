#!/usr/bin/env python3
# MIT License
#
# Copyright 2021 Broad Institute
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
Prepare and validate data by removing genes and covariates that would cause errors when running tensorQTL.
Optionally converts genotypes from BED to parquet format to reduce memory usage when reading in the genotypes.
"""

import argparse
from types import MethodType

import pandas as pd
import pyarrow as pa
import pyarrow.csv as csv
import pyarrow.parquet as pq

import misc_utils.pandas_utils as pd_utils
import dropseq_util.log_util as my_logging


def main() -> None:
    """
    Prepare and validate data by removing genes and covariates that would cause errors when running TensorQTL.
    Optionally converts genotypes from BED to parquet format to reduce memory usage when reading in the genotypes.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument('--genotype_bed', help='Genotypes BED file in tensorqtl format', required=False)
    parser.add_argument('--genotype_bim', help='Genotypes BIM file in PLINK format', required=False)
    parser.add_argument('--phenotypes', help='Phenotypes in BED format', required=False)
    parser.add_argument('--covariates', help='Covariates file, tab-delimited, covariates x samples', required=False)
    parser.add_argument('--genotype_out', help='Genotypes parquet output path', required=False)
    parser.add_argument('--phenotypes_out', help='Phenotypes output path', required=False)
    parser.add_argument('--covariates_out', help='Covariates output path', required=False)
    args = parser.parse_args()

    # If one of phenotypes, covariates, genotype_out, or phenotypes_out is provided, all must be provided
    if args.phenotypes or args.covariates or args.phenotypes_out or args.covariates_out:
        if not all([args.phenotypes, args.covariates, args.phenotypes_out, args.covariates_out]):
            parser.error('If one of --phenotypes, --covariates, --phenotypes_out, or --covariates_out is provided, '
                         'all must be provided')

        # If genotype_bim and genotype_bed are not provided, return a parser error
        if not args.genotype_bim and not args.genotype_bed:
            parser.error('Either --genotype_bim or --genotype_bed must be provided')

    # If no args specified, return a parser error
    if not args.phenotypes and not args.genotype_out:
        parser.error('Nothing to do, please specify at least one argument')

    if args.phenotypes:
        chrom_df = read_chrom_df(args.genotype_bim, args.genotype_bed)
        filter_tensorqtl_data(args.phenotypes, args.covariates, chrom_df, args.phenotypes_out, args.covariates_out)

    if args.genotype_bed and args.genotype_out:
        convert_bed_to_parquet(genotype_bed=args.genotype_bed, genotype_out=args.genotype_out)

    my_logging.log_message('Done')


def read_chrom_df(genotype_bim: str, genotype_bed: str) -> pd.DataFrame:
    """
    Read the chromosome dataframe from the genotype BIM file or BED file.
    """
    if genotype_bim:
        my_logging.log_message(f"Reading genotype chromosomes from '{genotype_bim}'")
        return pd.read_csv(genotype_bim, sep='\t', dtype=str, header=None, usecols=[0], names=['chrom'])
    else:
        my_logging.log_message(f"Reading genotype chromosomes from '{genotype_bed}'")
        return pd.read_csv(genotype_bed, sep='\t', dtype=str, header=0, usecols=[0], names=['chrom'])


def filter_tensorqtl_data(phenotypes: str,
                          covariates: str,
                          chrom_df: pd.DataFrame,
                          phenotypes_out: str,
                          covariates_out: str) -> None:
    my_logging.log_message(f"Reading phenotypes from '{phenotypes}'")
    phenotypes_df = pd.read_csv(phenotypes, sep='\t', dtype=str)
    my_logging.log_message(f"Reading covariates from '{covariates}'")
    covariates_df = pd.read_csv(covariates, sep='\t', comment='#', dtype=str, index_col=0, keep_default_na=False)

    covariates_cols = list(covariates_df.columns)
    covariates_rows = list(covariates_df.index)
    donors = list(phenotypes_df.columns[4:])

    # Set up the dataframe as samples x covariates for further processing
    if set(covariates_cols) & set(donors):
        # transpose the covariates for now, un-transpose later
        covariates_df = covariates_df.T
    elif set(covariates_rows) & set(donors):
        # rows are samples, so we're good for now, but will transpose the results later
        my_logging.log_message('Transposing covariates')
        pass
    else:
        raise Exception(
            f'Donors in {covariates} do not match donors in {phenotypes}:\n'
            f'Phenototypes donors: {donors}\n'
            f'Covariates rows: {covariates_rows}\n'
            f'Covariates columns: {covariates_cols}'
        )

    my_logging.log_message('Checking phenotypes and covariates')
    phenotypes_subset_df = phenotypes_df
    covariates_subset_df = covariates_df

    # Remove genes that do not have genotypes.
    phenotypes_pos_subset = phenotypes_subset_df.iloc[:, 0].isin(chrom_df['chrom'].unique())
    num_genes_no_genotype = phenotypes_subset_df.shape[0] - phenotypes_pos_subset.sum()
    if num_genes_no_genotype > 0:
        my_logging.log_message(f'Filtering out {num_genes_no_genotype} genes that do not overlap genotype chromosomes')
        phenotypes_subset_df = phenotypes_subset_df.loc[phenotypes_pos_subset]

    # List rows in phenotypes_subset_df where all values are the same as the first column's value
    phenotype_constant_rows = phenotypes_subset_df.iloc[:, 4:].eq(phenotypes_subset_df.iloc[:, 4], axis=0).all(axis=1)
    num_genes_const_expr = phenotype_constant_rows.sum()
    if num_genes_const_expr > 0:
        my_logging.log_message(f'Filtering out {num_genes_const_expr} genes with constant expression')
        phenotypes_subset_df = phenotypes_subset_df.loc[~phenotype_constant_rows]

    # Remove covariates for donors that were filtered out.
    num_covariates_no_donor = len(set(covariates_subset_df.index) - set(donors))
    if num_covariates_no_donor > 0:
        my_logging.log_message(f'Filtering out {num_covariates_no_donor} donors from the covariates')
        covariates_subset_df = covariates_subset_df.loc[donors]

    # Verify that covariates are all numeric.
    try:
        covariates_subset_df.astype(float)
    except ValueError as e:
        raise Exception(f'Covariates in {covariates} must be numeric') from e

    # If phenotypes_out ends in ".parquet", write as parquet, otherwise write as TSV.
    if phenotypes_out.endswith('.parquet'):
        my_logging.log_message(f"Writing phenotypes to '{phenotypes_out}'")
        # convert start and end columns to int32. convert donor columns to float64.
        first_column_types = {
            phenotypes_subset_df.columns[0]: 'string',
            phenotypes_subset_df.columns[1]: 'int32',
            phenotypes_subset_df.columns[2]: 'int32',
            phenotypes_subset_df.columns[3]: 'string',
        }
        donor_column_types = {
            donor: 'float64' for donor in donors
        }
        column_types = {**first_column_types, **donor_column_types}
        phenotypes_numeric_df = phenotypes_subset_df.astype(column_types)
        phenotypes_numeric_df.to_parquet(phenotypes_out, index=False)
    else:
        my_logging.log_message(f"Writing phenotypes to '{phenotypes_out}'")
        pd_utils.to_tsv(phenotypes_subset_df, phenotypes_out, index=False)
    # Write covariates output as covariates x samples
    pd_utils.to_tsv(covariates_subset_df.T, covariates_out, index=True)


def convert_bed_to_parquet(genotype_bed: str, genotype_out: str) -> None:
    """
    Convert a genotypes BED file to a binary parquet file.
    """

    # Read the first row to get the names of the columns using pandas.
    my_logging.log_message(f"Reading genotype data from '{genotype_bed}'")
    header_line = pd.read_csv(genotype_bed, sep='\t', dtype=str, nrows=0).columns

    first_fields = [
        pa.field(header_line[0], pa.string()),
        pa.field(header_line[1], pa.int32()),
        pa.field(header_line[2], pa.int32()),
        pa.field(header_line[3], pa.string()),
    ]
    donor_fields = [
        pa.field(header, pa.int8()) for header in header_line[4:]
    ]
    table_fields = first_fields + donor_fields

    read_options = csv.ReadOptions(block_size=2 ** 20)
    parse_options = csv.ParseOptions(delimiter='\t')
    convert_options = csv.ConvertOptions(column_types=pa.schema(table_fields))

    my_logging.log_message(f"Writing genotype data to '{genotype_out}'")
    table = csv.read_csv(
        input_file=genotype_bed,
        read_options=read_options,
        parse_options=parse_options,
        convert_options=convert_options,
    )
    pq.write_table(table=table, where=genotype_out)


if __name__ == '__main__':
    main()
