# MIT License
#
# Copyright 2026 Broad Institute
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
Optional QC step after scRNA-seq aggregation: a donor QC table and PDFs of library, donor and cell type QC.

Outputs are written to --output-dir, named <analysis-id>.<report>:
  library_metrics.pdf, filtered_library_metrics.pdf   always
  donor_qc.txt, donor_summary.pdf, donor_nuclei_counts.pdf   if the donors expected in the experiment are known
  cell_type_counts.pdf, cell_type_median_UMIs.pdf, cell_type_proportions.pdf   if cells have cell types
  donor_cell_type_distributions.pdf   if there are both

The expected donors are the donors in the --donor-file files, which are compared by content and must list the same
donors.  If there are none, or they differ, which is logged, the donor outputs are skipped.  Cells with no cell type are left out of the cell type pages.

The condition, an optional property of each library, comes from --conditions-file (a header and two columns, library
and condition) or from the cell metadata column --condition-column, which must have one value per library.

--core-cell-types is a headerless list of the cell types that are consistently present in the tissue, whatever the
dissection.  donor_qc.txt gets the median UMIs per nucleus of each donor for each of them present, and pooled.

Skipped outputs are not an error.  Invalid input exits with a non-zero status.
"""

import argparse
import os
import sys

import pandas as pd

from dropseq.aggregation import logger, add_log_argument, dctLogLevel
from dropseq.aggregation.plots import qc_tables
from dropseq.aggregation.plots.cell_type_reports import (cell_type_counts_report, cell_type_median_umis_report,
                                                         cell_type_proportions_report,
                                                         donor_cell_type_distributions_report)
from dropseq.aggregation.plots.donor_reports import donor_nuclei_counts_report, donor_summary_report
from dropseq.aggregation.plots.library_reports import filtered_library_metrics_report, library_metrics_report
from dropseq.aggregation.plots.qc_tables import CELL_TYPE, DONOR, LIBRARY, UMIS

class InvalidInput(Exception):
    pass


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    parser.add_argument("--cell-metadata", required=True, help="Tab-separated cell metadata of the aggregation.")
    parser.add_argument("--library-metrics", required=True, help="Tab-separated library metrics of the aggregation.")
    parser.add_argument("--analysis-id", required=True, help="Title of the pages and prefix of the output files.")
    parser.add_argument("--output-dir", default=".", help="Directory for the outputs. (default: %(default)s)")
    parser.add_argument("--donor-file", nargs="+", metavar="FILE",
                        help="Headerless file with the donors expected in the experiment, one per line, local or "
                             "gs://.  Give the donor file of each library; they must list the same donors.  Without "
                             "it there are no donor outputs.")
    parser.add_argument("--cell-type-column",
                        help="Cell metadata column with the cell type, for example "
                             "HMBA_Human_WB_v0.5_neighborhood_name.  Without it there are no cell type outputs.")
    conditions = parser.add_mutually_exclusive_group()
    conditions.add_argument("--conditions-file", help="Tab-separated file with a header and two columns, library "
                                                      "and condition.")
    conditions.add_argument("--condition-column", help="Cell metadata column with the condition of the library.")
    parser.add_argument("--core-cell-types", help="Headerless file with one core cell type per line.")
    parser.add_argument("--umi-column", default="num_retained_transcripts",
                        help="Cell metadata column with the number of UMIs. (default: %(default)s)")
    parser.add_argument("--donor-column", default="donor", help="Cell metadata column with the donor. "
                                                                 "(default: %(default)s)")
    parser.add_argument("--library-column", default="PREFIX", help="Cell metadata column with the library. "
                                                                    "(default: %(default)s)")
    parser.add_argument("--library-metrics-library-column", default="library",
                        help="Library metrics column with the library. (default: %(default)s)")
    return parser.parse_args(args)


def _require(table, columns, source):
    missing = [c for c in columns if c not in table.columns]
    if missing:
        raise InvalidInput(f"No {missing} column(s) in {source}")


def run(options):
    os.makedirs(options.output_dir, exist_ok=True)
    out = options.output_dir
    analysis_id = options.analysis_id

    string_columns = {c: str for c in [options.library_column, options.donor_column, options.cell_type_column,
                                       options.condition_column] if c}
    cell_metadata = pd.read_csv(options.cell_metadata, sep="\t", dtype=string_columns)
    _require(cell_metadata, [options.library_column, options.umi_column, options.condition_column or
                             options.library_column], options.cell_metadata)
    library_metrics = pd.read_csv(options.library_metrics, sep="\t",
                                  dtype={options.library_metrics_library_column: str})
    _require(library_metrics, [options.library_metrics_library_column], options.library_metrics)
    library_metrics = library_metrics.rename(columns={options.library_metrics_library_column: LIBRARY})
    if cell_metadata[options.library_column].isna().any():
        raise InvalidInput(f"Some cells have no {options.library_column} in {options.cell_metadata}")

    # gates
    has_donors = options.donor_column in cell_metadata.columns and cell_metadata[options.donor_column].notna().any()
    if not has_donors:
        logger.info(f"No cell has a {options.donor_column}; skipping donor QC")
    has_cell_types = False
    if options.cell_type_column is None:
        logger.info("No --cell-type-column; skipping cell type QC")
    elif options.cell_type_column not in cell_metadata.columns:
        logger.info(f"No {options.cell_type_column} column in {options.cell_metadata}; skipping cell type QC")
    elif cell_metadata[options.cell_type_column].isna().all():
        logger.info(f"No cell has a {options.cell_type_column}; skipping cell type QC")
    else:
        has_cell_types = True
        unlabelled = int(cell_metadata[options.cell_type_column].isna().sum())
        if unlabelled:
            logger.info(f"{unlabelled} of {len(cell_metadata)} cells have no cell type; they are left out of the "
                        "cell type pages")

    # inputs
    try:
        conditions = qc_tables.resolve_conditions(
            cell_metadata, options.library_column, options.conditions_file, options.condition_column,
            libraries=library_metrics[LIBRARY])
        expected_donors = qc_tables.expected_donors_from_files(options.donor_file)
        core_cell_types = qc_tables.read_core_cell_types(options.core_cell_types) if options.core_cell_types else []
    except (ValueError, OSError) as e:
        raise InvalidInput(str(e))
    do_donors = has_donors and expected_donors is not None
    if core_cell_types and not (do_donors and has_cell_types):
        logger.info("Ignoring --core-cell-types, which need both donor QC and cell types")
        core_cell_types = []

    cells = qc_tables.normalize_cell_metadata(
        cell_metadata, options.library_column, options.donor_column,
        options.cell_type_column if has_cell_types else None, options.umi_column,
        conditions)
    typed_cells = cells.dropna(subset=[CELL_TYPE]) if has_cell_types else None

    # library level
    cell_summary = (qc_tables.library_cell_summary(cells, CELL_TYPE, LIBRARY, UMIS) if has_cell_types else None)
    library_metrics_report(library_metrics, conditions, analysis_id, out, cell_summary)
    filtered_library_metrics_report(qc_tables.library_cell_summary(cells, None, LIBRARY, UMIS), conditions,
                                    analysis_id, out)

    # donor level
    donor_cell_type = None
    if do_donors:
        donor_qc = qc_tables.donor_qc_table(cells, expected_donors, core_cell_types,
                                            CELL_TYPE if has_cell_types else None, DONOR, UMIS)
        donor_qc.assign(expected=donor_qc["expected"].map({True: "true", False: "false"})).to_csv(
            os.path.join(out, f"{analysis_id}.donor_qc.txt"), sep="\t", index=False, na_rep="NA",
            float_format="%.15g")
        if has_cell_types:
            donor_cell_type = qc_tables.complete_donor_cell_type(typed_cells, CELL_TYPE, expected_donors, DONOR, UMIS)
        donor_summary_report(donor_qc, donor_cell_type, analysis_id, out)
        equitability = None
        if "donor_equitability" in library_metrics.columns:
            equitability = dict(zip(library_metrics[LIBRARY].astype(str), library_metrics["donor_equitability"]))
            equitability = {k: v for k, v in equitability.items() if pd.notna(v)}
        donor_nuclei_counts_report(cells, donor_qc, analysis_id, out, equitability)

    # cell type level
    if has_cell_types:
        cell_type_counts_report(typed_cells, donor_cell_type, analysis_id, out)
        cell_type_median_umis_report(typed_cells, donor_cell_type, analysis_id, out)
        cell_type_proportions_report(typed_cells, donor_cell_type, analysis_id, out)
        if donor_cell_type is not None:
            donor_cell_type_distributions_report(typed_cells, donor_cell_type, analysis_id, out)


def main(args=None):
    options = parse_args(args)
    logger.setLevel(dctLogLevel[options.log_level])
    try:
        run(options)
    except InvalidInput as e:
        logger.error(str(e))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
