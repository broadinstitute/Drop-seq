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
Convert optimus h5ad to standard Drop-seq output files.
"""
import argparse
import anndata as ad
import pandas as pd
from scipy.io import mmwrite

try:
    from . import cli
    from . import io_utils
except ImportError:
    import cli
    import io_utils


def add_subparser(subparsers):
    parser = subparsers.add_parser("optimus_h5ad_to_dropseq", description=__doc__)
    parser.add_argument("--input", "-i", required=True, help="Input h5ad file.")
    parser.add_argument("--min-transcripts", type=int, default=20,
                        help="Minimum number of transcripts to keep a CBC.  (default: %(default)s)")
    parser.add_argument("--h5ad",
                        help="Output same format as input, but with CBCs < min-transcripts and duplicated genes "
                             "removed.")
    parser.add_argument("--mtx", help="Output mtx. If .gz extension, will be gzipped.")
    parser.add_argument("--barcodes", help="Output barcodes. If .gz extension, will be gzipped.")
    parser.add_argument("--features", help="Output features. If .gz extension, will be gzipped.")
    parser.add_argument("--dge", help="Output tab-separated text file, with rows=genes and columns=cells. "
                                      "If .gz extension, will be gzipped.")
    parser.add_argument("--summary", type=argparse.FileType(mode="w"),
                        help="Output a tab-separated DGE summary text file with columns CELL_BARCODE, NUM_GENIC_READS, "
                             "NUM_TRANSCRIPTS, NUM_GENES.")
    parser.add_argument("--reads-per-cell",
                        help="Output a tab-separated text file with columns num_reads, cell_names.")
    parser.add_argument("--read-quality-metrics", type=argparse.FileType(mode="w"),
                        help="Output a tab-separated ReadQualityMetrics text file.")
    parser.add_argument("--cell-selection-report", help="Output a table of per-cell-barcode metrics.")


def main(options):
    cli.logger.info(f'loading full adata {options.input}')
    adata = ad.read_h5ad(options.input)

    duplicated_gene_name = adata.var['gene_names'].duplicated(keep=False)
    num_duplicated_genes = duplicated_gene_name.sum()
    if num_duplicated_genes > 0:
        cli.logger.info(f'Removing {num_duplicated_genes} duplicated gene names')
        adata = adata[:, ~duplicated_gene_name]

    # Convert the counts matrix to integers, and transpose where rows are genes and columns are cell barcodes
    matrix = adata.X.astype(int).T
    # Count the number of transcripts per cell barcode
    num_transcripts = matrix.sum(axis=0).A1

    cli.logger.info(f'subsetting to barcodes with at least {options.min_transcripts} transcripts')
    adata = adata[num_transcripts >= options.min_transcripts, :]

    # Convert the counts matrix to integers, and transpose where rows are genes and columns are cell barcodes
    matrix = adata.X.astype(int).T
    # Count the number of transcripts per cell barcode
    num_transcripts = matrix.sum(axis=0).A1

    # add additional columns with names that are expected by the downstream tools
    obs = adata.obs.copy()
    total_reads = obs['n_reads']
    mapped_reads = obs['reads_mapped_uniquely']
    obs['NUM_GENES'] = obs['n_genes']
    obs['NUM_GENIC_READS'] = obs['reads_mapped_exonic'] + obs['reads_mapped_exonic_as'] + \
                             obs['reads_mapped_intronic'] + obs['reads_mapped_intronic_as']
    obs['NUM_TRANSCRIPTS'] = num_transcripts
    obs['num_transcripts'] = num_transcripts
    obs['num_reads'] = mapped_reads
    obs['totalReads'] = total_reads
    obs['mappedReads'] = mapped_reads
    obs['hqMappedReads'] = mapped_reads
    obs['hqMappedReadsNoPCRDupes'] = mapped_reads
    obs['pct_coding'] = (obs['reads_mapped_exonic'] + obs['reads_mapped_exonic_as']) / mapped_reads
    obs['pct_intronic'] = (obs['reads_mapped_intronic'] + obs['reads_mapped_intronic_as']) / mapped_reads
    obs['pct_intergenic'] = obs['reads_mapped_intergenic'] / mapped_reads
    obs['pct_mt'] = obs['reads_mapped_mitochondrial'] / mapped_reads
    obs['pct_genic'] = obs['pct_coding'] + obs['pct_intronic']
    obs['pct_ribosomal'] = 0
    # NOTE: Unlike Drop-seq, Optimus does not generate separate UTR counts from coding counts.
    # Thus Optimus[reads_mapped_exonic + reads_mapped_exonic_as] == Drop-seq[coding + UTR]
    # It may be more correct to generate pct_utr in the cell_selection_report as NA values,
    # but that might require changes to the downstream tools.
    obs['pct_utr'] = 0

    obs['pct_coding'] = obs['pct_coding'].round(4)
    obs['pct_intronic'] = obs['pct_intronic'].round(4)
    obs['pct_intergenic'] = obs['pct_intergenic'].round(4)
    obs['pct_mt'] = obs['pct_mt'].round(4)
    obs['pct_genic'] = obs['pct_genic'].round(4)
    obs['pct_ribosomal'] = obs['pct_ribosomal'].round(4)
    obs['pct_utr'] = obs['pct_utr'].round(4)

    if options.h5ad is not None:
        cli.logger.info(f'generating h5ad')
        adata.write(options.h5ad)

    if options.mtx is not None:
        cli.logger.info('generating mtx')
        with io_utils.open_maybe_gz(options.mtx, 'wb') as fOut:
            mmwrite(fOut, matrix)
    if options.barcodes is not None:
        cli.logger.info(f'generating barcodes')
        with io_utils.open_maybe_gz(options.barcodes, 'wt') as fOut:
            adata.obs_names.to_series().to_csv(fOut, header=False, index=False)
    if options.features is not None:
        cli.logger.info(f'generating features')
        features_df = pd.DataFrame(adata.var_names.to_series(), columns=['gene_id'])
        features_df['gene_name'] = features_df['gene_id']
        features_df['feature_type'] = 'Gene Expression'
        with io_utils.open_maybe_gz(options.features, 'wt') as fOut:
            features_df.to_csv(fOut, sep='\t', header=False, index=False)

    if options.dge is not None:
        cli.logger.info(f'generating dge')
        dge = pd.DataFrame.sparse.from_spmatrix(matrix)
        dge.columns = adata.obs_names
        dge.index = adata.var_names
        dge.index.name = 'GENE'
        with io_utils.open_maybe_gz(options.dge, 'wt') as fOut:
            fOut.write('#DGE\tVERSION:1.1\tEXPRESSION_FORMAT:raw\n')
            dge.to_csv(fOut, sep='\t')

    if options.summary is not None:
        cli.logger.info(f'generating summary')
        dge_summary = obs[['NUM_GENIC_READS', 'NUM_TRANSCRIPTS', 'NUM_GENES']]
        dge_summary.index.name = 'CELL_BARCODE'
        dge_summary = dge_summary.sort_values(by='NUM_GENIC_READS', ascending=False)
        options.summary.write('## METRICS CLASS\torg.broadinstitute.dropseqrna.barnyard.DigitalExpression$DESummary\n')
        dge_summary.to_csv(options.summary, sep='\t')

    if options.reads_per_cell is not None:
        with io_utils.open_maybe_gz(options.reads_per_cell, 'wt') as fOut:
            cli.logger.info('generating reads per cell')
            reads_per_cell = obs[['num_reads', 'cell_names']]
            reads_per_cell = reads_per_cell.sort_values(by='num_reads', ascending=False)
            reads_per_cell.to_csv(fOut, sep='\t', header=False, index=False)

    if options.read_quality_metrics is not None:
        cli.logger.info('generating read quality metrics')
        read_qualities = obs[['totalReads', 'mappedReads', 'hqMappedReads', 'hqMappedReadsNoPCRDupes']]
        read_quality_metrics = pd.DataFrame(read_qualities.sum()).T
        read_quality_metrics.insert(0, 'aggregate', 'all')
        options.read_quality_metrics.write(
            '## METRICS CLASS\torg.broadinstitute.dropseqrna.metrics.ReadQualityMetrics\n')
        read_quality_metrics.to_csv(options.read_quality_metrics, sep='\t', index=False)

    if options.cell_selection_report is not None:
        cli.logger.info('generating cell selection report')
        cell_selection_report = obs[
            ['num_transcripts', 'num_reads', 'pct_ribosomal', 'pct_coding', 'pct_intronic', 'pct_intergenic', 'pct_utr',
             'pct_genic', 'pct_mt']]
        cell_selection_report.index.name = 'cell_barcode'
        cell_selection_report = cell_selection_report.sort_values(by='num_transcripts', ascending=False)
        cell_selection_report.to_csv(options.cell_selection_report, sep='\t')

    cli.logger.info('Done')
    return 0
