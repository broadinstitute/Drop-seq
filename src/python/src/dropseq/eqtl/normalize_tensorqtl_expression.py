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
Normalize tensorQTL expression per gene with edgeR CPM then per donor with an inverse normal transformation.
"""

# Based on a combination of:
#  - https://github.com/broadinstitute/pyqtl/blob/v0.1.8/qtl/norm.py
#  - https://github.com/broadinstitute/eqtl_pipeline_terra/blob/cc7eca7/dockerfiles/preprocess/src/python/normalize.py
#
# Changes:
#  - Use math.log(dbl, 2) instead of numpy.log2(dbl) so that singularity outputs match non-singularity outputs
#  - Input is also in tensorQTL's phenotypes BED format

import argparse
import math
import sys
import warnings
from typing import Optional, List

import numpy as np
import pandas as pd
import scipy.stats as stats



def consistent_log2(dbl):
    """
    Use the slower and possibly loss of precision math.log(dbl, 2) instead of np.log2 unless math.log cannot process the
    value.

    This handles cases where on some machines these functions may produce different outputs for certain values of dbl:
     - np.log2(dbl)
     - math.log2(dbl)
     - math.log(dbl, 2)

    One example input: dbl = 5.96392877939081911

    Depending on the environment, print(f'{log2(d):.20f}') is either 2.57626302921203498286 or 2.57626302921203542695
    for functions np.log2 and math.log2.

    np.log2 and math.log2 may produce different results even on the same host. When running on uger nodes,
    sv01, in docker on one's x86 macbook, or a GCE instance they may produce different outputs on the same host
    depending on if they are executed inside or outside a docker/singularity container.

    Meanwhile, math.log(dbl, 2) is consistent as far as I can tell. For the above example, it always produces
    2.57626302921203542695 on all hosts / environments, so far.

    May or may not have something to do with issues discussed in these links:
     - https://stackoverflow.com/questions/17702065/python-numpy-log2-vs-matlab#answer-17702094
     - https://github.com/numpy/numpy/issues/4787
     - https://github.com/numpy/numpy/issues/13836
     - https://github.com/python/cpython/issues/47974
    """
    try:
        return math.log(dbl, 2)
    except ValueError:
        return np.log2(dbl)


# Modified from https://github.com/broadinstitute/pyqtl/blob/v0.1.8/qtl/norm.py#L104-L165
# noinspection PyPep8Naming
def edger_calcnormfactors(counts_df, ref=None, logratio_trim=0.3,
                          sum_trim=0.05, acutoff=-1e10, verbose=False):
    """
    Calculate TMM (Trimmed Mean of M values) normalization.
    Reproduces edgeR::calcNormFactors.default
    Scaling factors for the library sizes that minimize
    the log-fold changes between the samples for most genes.
    Effective library size: TMM scaling factor * library size
    References:
     [1] Robinson & Oshlack, 2010
     [2] R functions:
          edgeR::calcNormFactors.default
          edgeR:::.calcFactorWeighted
          edgeR:::.calcFactorQuantile
    """

    # discard genes with all-zero counts
    Y = counts_df.values.copy()
    allzero = np.sum(Y > 0, axis=1) == 0
    if np.any(allzero):
        Y = Y[~allzero, :]

    # select reference sample
    if ref is None:  # reference sample index
        f75 = np.percentile(Y/np.sum(Y, axis=0), 75, axis=0)
        ref = np.argmin(np.abs(f75-np.mean(f75)))
        if verbose:
            print('Reference sample index: '+str(ref))

    N = np.sum(Y, axis=0)  # total reads in each library

    # (Mostly) use a vectorized math.log2 instead of np.log2
    vec_consistent_log2 = np.vectorize(consistent_log2)

    # with np.errstate(divide='ignore'):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # log fold change; Mg in [1]
        logR = vec_consistent_log2((Y/N).T / (Y[:, ref]/N[ref])).T
        # average log relative expression; Ag in [1]
        absE = 0.5*(vec_consistent_log2(Y/N).T + vec_consistent_log2(Y[:, ref]/N[ref])).T
        v = (N-Y)/N/Y
        v = (v.T + v[:, ref]).T  # w in [1]

    ns = Y.shape[1]
    tmm = np.zeros(ns)
    for i in range(ns):
        fin = np.isfinite(logR[:, i]) & np.isfinite(absE[:, i]) & (absE[:, i] > acutoff)
        n = np.sum(fin)

        loL = np.floor(n*logratio_trim)+1
        hiL = n + 1 - loL
        loS = np.floor(n*sum_trim)+1
        hiS = n + 1 - loS
        rankR = stats.rankdata(logR[fin, i])
        rankE = stats.rankdata(absE[fin, i])
        keep = (rankR >= loL) & (rankR <= hiL) & (rankE >= loS) & (rankE <= hiS)
        # in [1], w erroneously defined as 1/v ?
        tmm[i] = 2**(np.nansum(logR[fin, i][keep]/v[fin, i][keep]) / np.nansum(1/v[fin, i][keep]))

    tmm = tmm / np.exp(np.mean(np.log(tmm)))
    return tmm


# Modified from https://github.com/broadinstitute/pyqtl/blob/v0.1.8/qtl/norm.py#L186-L197
def edger_cpm(counts_df, tmm=None, normalized_lib_sizes=True):
    """
    Return edgeR normalized/rescaled CPM (counts per million)
    Reproduces edgeR::cpm.DGEList
    """

    lib_size = counts_df.sum(axis=0)
    if normalized_lib_sizes:
        if tmm is None:
            tmm = edger_calcnormfactors(counts_df)
        lib_size = lib_size * tmm
    return counts_df / lib_size * 1e6


# Copied from https://github.com/broadinstitute/pyqtl/blob/v0.1.8/qtl/norm.py#L57-L67
# noinspection PyPep8Naming
def inverse_normal_transform(M):
    """Transform rows to a standard normal distribution"""
    if isinstance(M, pd.Series):
        r = stats.mstats.rankdata(M)
        return pd.Series(stats.norm.ppf(r/(M.shape[0]+1)), index=M.index, name=M.name)
    else:
        R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
        Q = stats.norm.ppf(R/(M.shape[1]+1))
        if isinstance(M, pd.DataFrame):
            Q = pd.DataFrame(Q, index=M.index, columns=M.columns)
        return Q


# Modified from
# https://github.com/broadinstitute/eqtl_pipeline_terra/blob/cc7eca7/dockerfiles/preprocess/src/python/normalize.py
def main(args: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        dest='phenotype_bed',
        type=str,
        help='Phenotypes in BED format. '
             'For more information see https://github.com/broadinstitute/tensorqtl/tree/v1.0.7#input-formats',
    )
    parser.add_argument(
        dest='output_prefix',
        type=str,
        help='prefix for output files',
    )
    options = parser.parse_args(args)

    # read in genes x donors count matrix
    phenotype_df = pd.read_csv(options.phenotype_bed, sep='\t', index_col=None)

    phenotype_mapping = {
        phenotype_df.columns[0]: '#chr',
        phenotype_df.columns[1]: 'start',
        phenotype_df.columns[2]: 'end',
        phenotype_df.columns[3]: 'pid',
    }
    phenotype_df = phenotype_df.rename(columns=phenotype_mapping)

    # sort [chr1, chr10,..chr2, chr20,.., chr3,..chr9]
    phenotype_df = phenotype_df.sort_values(['#chr', 'start']).set_index('pid', drop=False)

    # edgeR CPM normalization
    cpm_df = edger_cpm(phenotype_df.iloc[:, 4:])

    out_df = phenotype_df.iloc[:, :4].join(cpm_df)
    out_df = out_df.sort_values(['#chr', 'start'])

    out_file = f'{options.output_prefix}.TPM_expression'
    out_df.to_csv(out_file + '.bed', sep='\t', index=False)

    # inverse normal transform
    int_df = inverse_normal_transform(cpm_df)

    out_df = phenotype_df.iloc[:, :4].join(int_df)
    out_df = out_df.sort_values(['#chr', 'start'])

    out_file = f'{options.output_prefix}.normalized_expression'
    out_df.to_csv(out_file + '.bed', sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
