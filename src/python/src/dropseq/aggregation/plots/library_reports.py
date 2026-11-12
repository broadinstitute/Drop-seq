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
Library-level QC reports: <id>.library_metrics.pdf and <id>.filtered_library_metrics.pdf.
"""

import os
from collections import namedtuple

from dropseq.aggregation import logger
from dropseq.aggregation.plots import plot_primitives as pp
from dropseq.aggregation.plots.qc_tables import CONDITION, LIBRARY

LibraryMetric = namedtuple("LibraryMetric", ["column", "label", "y_max"], defaults=[None])

# The pages of the library metrics report, in order.  A column that is missing or all NA is skipped.
LIBRARY_METRICS = [
    LibraryMetric("total_reads_M", "Total reads (M)"),
    LibraryMetric("pct_hq_mapped_reads", "% uniquely mapped reads", 100),
    LibraryMetric("pct_chimeric_UMIs", "% chimeric UMIs"),
    LibraryMetric("n_STAMPs", "# nuclei SELECTED"),
    LibraryMetric("median_cell_UMIs_cbrb", "Median UMIs/selected nuclei (post-CBRB)"),
    LibraryMetric("pct_confident_doublets", "% donor doublets"),
    LibraryMetric("assignable_singlets", "# nuclei assigned to a donor"),
    LibraryMetric("reads_per_umi_donor", "Reads/UMI (in assignable singlets)"),
    LibraryMetric("totalUMIsM", "Total UMIs (millions, in assignable singlets)"),
    LibraryMetric("scDblFinder_pct_doublets", "% scDblFinder doublets"),
]

# Pages from the cell metadata, drawn only if there are cell types: columns of qc_tables.library_cell_summary.
CELL_TYPE_LIBRARY_METRICS = [
    LibraryMetric("num_labelled_nuclei", "# nuclei with a cell type"),
    LibraryMetric("median_umis_per_labelled_nucleus", "Median UMIs/nucleus with a cell type"),
]

# The library metrics pages are narrower than the others, as in the R report.
LIBRARY_PAGE_WIDTH = 10
LIBRARY_SUBTITLE = lambda metric: f"{metric.label} per library. Dashed line is median."

NO_CONDITION = "same conditions"
ALL_LIBRARIES = "All libraries"


def with_conditions(table, conditions, library_column=LIBRARY, no_condition=NO_CONDITION):
    """
    :return: copy of table (which has a column library_column) with a condition column, which is no_condition for
    every library if conditions is None
    """
    table = table.rename(columns={library_column: LIBRARY}) if library_column != LIBRARY else table.copy()
    table[LIBRARY] = table[LIBRARY].astype(str)
    if conditions is None:
        table[CONDITION] = no_condition
    else:
        table[CONDITION] = table[LIBRARY].map(dict(zip(conditions[LIBRARY], conditions[CONDITION]))).fillna("unknown")
    return table


def _metric_pages(table, metrics, analysis_id, subtitle, width, median):
    """
    The page of each metric that is available: a bar per library (coloured as in R, with ggplot's default hue
    palette) and a box by condition.

    :param subtitle: function from a LibraryMetric to the subtitle of its page
    """
    pages = []
    for metric in metrics:
        if metric.column not in table.columns or table[metric.column].isna().all():
            logger.info(f"Skipping library metric {metric.column}: missing or all NA")
            continue
        pages.append(lambda m=metric: pp.bar_and_condition_box(
            table, m.column, m.label, None, analysis_id, subtitle(m), y_max=m.y_max, width=width, median=median))
    return pages


def library_metrics_report(library_metrics, conditions, analysis_id, out_dir, cell_summary=None,
                           library_column="library"):
    """
    Write <analysis_id>.library_metrics.pdf, a page per metric of LIBRARY_METRICS that is available.

    :param library_metrics: DataFrame with one row per library
    :param conditions: library -> condition table, or None
    :param cell_summary: qc_tables.library_cell_summary of the cells with their cell type, or None if there are no
    cell types, which adds the CELL_TYPE_LIBRARY_METRICS pages
    :return: the path, or None if no metric was available
    """
    table = library_metrics.copy()
    if "totalUMIs" in table.columns:
        table["totalUMIsM"] = (table["totalUMIs"] / 1e6).round(2)
    table = with_conditions(table, conditions, library_column)
    pages = _metric_pages(table, LIBRARY_METRICS, analysis_id, LIBRARY_SUBTITLE, LIBRARY_PAGE_WIDTH, True)
    if cell_summary is not None:
        summary = with_conditions(cell_summary, conditions)
        pages += _metric_pages(summary, CELL_TYPE_LIBRARY_METRICS, analysis_id, LIBRARY_SUBTITLE,
                               LIBRARY_PAGE_WIDTH, True)
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.library_metrics.pdf"))


def filtered_library_metrics_report(library_summary, conditions, analysis_id, out_dir):
    """
    Write <analysis_id>.filtered_library_metrics.pdf: the number of nuclei and their median UMIs per library, after
    aggregation filtering.

    :param library_summary: qc_tables.library_cell_summary
    """
    table = with_conditions(library_summary, conditions, no_condition=ALL_LIBRARIES)
    metrics = [LibraryMetric("num_nuclei", "# nuclei"),
               LibraryMetric("median_umis_per_nucleus", "Median UMIs/nucleus")]
    subtitles = {"num_nuclei": "Nuclei counts per library.",
                 "median_umis_per_nucleus": "Median UMIs/nucleus per library."}
    pages = _metric_pages(table, metrics, analysis_id, lambda m: subtitles[m.column], pp.PAGE_SIZE[0], False)
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.filtered_library_metrics.pdf"))
