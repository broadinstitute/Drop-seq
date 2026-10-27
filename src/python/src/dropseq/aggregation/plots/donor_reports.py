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
Donor QC reports: <id>.donor_summary.pdf and <id>.donor_nuclei_counts.pdf.
"""

import os

from dropseq.aggregation.plots import plot_primitives as pp
from dropseq.aggregation.plots.cell_type_reports import (COUNTS, MEDIAN_UMIS, PROPORTIONS,
                                                         donor_cell_type_distribution_pages,
                                                         donor_cell_type_pages)
from dropseq.aggregation.plots.qc_tables import (CONDITION, CORE_POOLED_COLUMN, CORE_TYPE_COLUMN_PREFIX, DONOR,
                                                 LIBRARY, donor_library_counts)

DONOR_LABEL = "Donor"


def _donor_bar_page(donor_qc, column, ylabel, subtitle, colors, analysis_id):
    unexpected = list(donor_qc.loc[~donor_qc["expected"], DONOR])
    return lambda: pp.sorted_bar_with_median(donor_qc, DONOR, column, ylabel, colors, analysis_id, subtitle,
                                             red=unexpected)


def donor_summary_report(donor_qc, donor_cell_type, analysis_id, out_dir):
    """
    Write <analysis_id>.donor_summary.pdf: bars of the donor QC table per donor and, if donor_cell_type (from
    qc_tables.complete_donor_cell_type) is given, the cell type composition and quality of each donor.

    :param donor_qc: qc_tables.donor_qc_table
    """
    colors = pp.donor_palette(donor_qc[DONOR])
    table = donor_qc.assign(percent_nuclei=(donor_qc["fraction_nuclei"] * 100).round(3),
                            total_umis_m=donor_qc["total_umis_m"].round(2),
                            median_umis_per_nucleus=donor_qc["median_umis_per_nucleus"].round(0))
    page = lambda column, ylabel, subtitle: _donor_bar_page(table, column, ylabel, subtitle, colors, analysis_id)
    pages = [
        page("num_nuclei", "# nuclei", "Total nuclei counts per donor."),
        page("percent_nuclei", "% of nuclei", "Percent of nuclei per donor."),
        page("total_umis_m", "Total # of UMIs (M)", "Total UMI counts per donor."),
        page("median_umis_per_nucleus", "Median UMIs / nucleus",
             "Median UMIs / nucleus per donor (rounded to nearest integer). Donors with at least one nucleus."),
    ]
    if CORE_POOLED_COLUMN in table.columns:
        pooled = table.assign(**{CORE_POOLED_COLUMN: table[CORE_POOLED_COLUMN].round(0)})
        pages.append(_donor_bar_page(pooled, CORE_POOLED_COLUMN, "Median UMIs / core cell type nucleus",
                                     "Median UMIs / nucleus of all core cell types combined, per donor "
                                     "(rounded to nearest integer).", colors, analysis_id))
        for column in [c for c in table.columns if c.startswith(CORE_TYPE_COLUMN_PREFIX)]:
            core_type = column[len(CORE_TYPE_COLUMN_PREFIX):]
            per_type = table.assign(**{column: table[column].round(0)})
            pages.append(_donor_bar_page(per_type, column, f"Median UMIs / {core_type}",
                                         f"Median UMIs / {core_type} nucleus per donor (rounded to nearest integer). "
                                         f"NA if no {core_type} nucleus was detected.", colors, analysis_id))
    if donor_cell_type is not None:
        pages += donor_cell_type_distribution_pages(donor_cell_type, analysis_id)
        for metric in (COUNTS, PROPORTIONS, MEDIAN_UMIS):
            pages += donor_cell_type_pages(donor_cell_type, metric, colors, analysis_id)
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.donor_summary.pdf"))


def donor_nuclei_counts_report(cells, donor_qc, analysis_id, out_dir, donor_equitability=None):
    """
    Write <analysis_id>.donor_nuclei_counts.pdf: the nuclei of each donor in total, by condition and by library.

    :param cells: normalized cell metadata
    :param donor_qc: qc_tables.donor_qc_table
    :param donor_equitability: optional dict from library to its donor equitability, printed above the library's bar
    of the donor composition page
    """
    colors = pp.donor_palette(donor_qc[DONOR])
    has_conditions = CONDITION in cells.columns
    counts = donor_library_counts(cells)
    library_colors = pp.palette_for(counts[LIBRARY])
    donor_order = list(donor_qc.sort_values("num_nuclei", ascending=False, kind="stable")[DONOR])
    pages = [_donor_bar_page(donor_qc, "num_nuclei", "# nuclei", "Total nuclei counts per donor.", colors,
                             analysis_id)]
    if has_conditions:
        condition_colors = pp.palette_for(counts[CONDITION])
        pages.append(lambda: pp.stacked_bars(counts, DONOR, CONDITION, "num_nuclei", condition_colors, "# of nuclei",
                                             analysis_id, "Donor nuclei counts per condition.", DONOR_LABEL,
                                             normalize=False, label_min=0, x_order=donor_order))
    pages.append(lambda: pp.box_with_points(counts, DONOR, "num_nuclei", "# of nuclei", analysis_id,
                                            "Each point is a library.", DONOR_LABEL, color_by=LIBRARY,
                                            colors=library_colors, order=donor_order, legend_title="Library"))
    if has_conditions:
        pages.append(lambda: pp.box_with_points(counts, DONOR, "num_nuclei", "# of nuclei", analysis_id,
                                                "Each point is a library.", DONOR_LABEL, hue=CONDITION,
                                                hue_colors=condition_colors, order=donor_order,
                                                legend_title="Condition"))
    library_order = list(counts.groupby(LIBRARY)["num_nuclei"].sum().sort_values(ascending=False).index)
    pages.append(lambda: pp.box_with_points(counts, LIBRARY, "num_nuclei", "# of nuclei", analysis_id,
                                            "Each point is a donor.", "Library", color_by=DONOR, colors=colors,
                                            facet=CONDITION if has_conditions else None, order=library_order,
                                            legend_title=DONOR_LABEL))
    annotations = None
    subtitle = "Donor composition per library."
    if donor_equitability:
        annotations = {library: f"{value:.2f}" for library, value in donor_equitability.items()}
        subtitle += " Number above each bar is the donor equitability."
    pages.append(lambda: pp.stacked_bars(counts, LIBRARY, DONOR, "num_nuclei", colors, "% of nuclei", analysis_id,
                                         subtitle, "Library", facet=CONDITION if has_conditions else None,
                                         annotations=annotations))
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.donor_nuclei_counts.pdf"))
