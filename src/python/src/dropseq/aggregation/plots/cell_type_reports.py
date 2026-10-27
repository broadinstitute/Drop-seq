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
Cell type QC reports: <id>.cell_type_counts.pdf, <id>.cell_type_median_UMIs.pdf, <id>.cell_type_proportions.pdf and
<id>.donor_cell_type_distributions.pdf.

The cells are the normalized cell metadata of qc_tables.normalize_cell_metadata, which has CELL_TYPE and, if there are
conditions, CONDITION.  The donor_cell_type table is qc_tables.complete_donor_cell_type, or None if there is no donor
QC.
"""

import os
from collections import namedtuple

import pandas as pd

from dropseq.aggregation.plots import plot_primitives as pp
from dropseq.aggregation.plots.qc_tables import (CELL_TYPE, CONDITION, DONOR, LIBRARY, cell_type_table)

# column of cell_type_table, label of the value, factor to multiply it by, description of the pages
CellTypeMetric = namedtuple("CellTypeMetric", ["column", "label", "scale", "noun"])

COUNTS = CellTypeMetric("num_nuclei", "# nuclei", 1, "counts")
PROPORTIONS = CellTypeMetric("fraction_nuclei", "% of nuclei", 100, "proportions")
MEDIAN_UMIS = CellTypeMetric("median_umis_per_nucleus", "Median UMIs/nucleus", 1, "median UMIs/nucleus")

VALUE = "value"
CELL_TYPE_LABEL = "Cell type"


def _scaled(table, metric):
    table = table.copy()
    table[VALUE] = table[metric.column] * metric.scale
    if metric is MEDIAN_UMIS:
        table[VALUE] = table[VALUE].round(0)
    return table


def _with_condition(table, cells):
    """
    :return: table, which has a LIBRARY column, with the condition of each library, if the cells have one
    """
    if CONDITION in cells.columns:
        conditions = cells.drop_duplicates(LIBRARY).set_index(LIBRARY)[CONDITION]
        table = table.assign(**{CONDITION: table[LIBRARY].map(conditions)})
    return table


def donor_cell_type_pages(donor_cell_type, metric, donor_colors, analysis_id):
    """
    The pages of a metric per donor and cell type: a box plot over donors for each cell type, and a panel per cell
    type with a bar per donor.  Each is a function that returns the page.
    """
    table = _scaled(donor_cell_type, metric)
    order = list(table.dropna(subset=[VALUE]).groupby(CELL_TYPE)[VALUE].median().sort_values(ascending=False).index)
    noun = metric.noun
    return [
        lambda: pp.box_with_points(table, CELL_TYPE, VALUE, metric.label, analysis_id,
                                   f"Overall cell type {noun}. Each point is a donor.", CELL_TYPE_LABEL,
                                   color_by=DONOR, colors=donor_colors, order=order, legend_title="Donor"),
        lambda: pp.faceted_sorted_bars(table, CELL_TYPE, DONOR, VALUE, donor_colors, metric.label, analysis_id,
                                       f"Cell type {noun} per donor.", red_by="expected"),
    ]


def library_pages(cells, metric, analysis_id, order):
    """
    The pages of a metric per cell type of the libraries, and of the conditions if the cells have them.
    """
    by_library = _scaled(_with_condition(cell_type_table(cells, [LIBRARY]), cells), metric)
    library_colors = pp.palette_for(by_library[LIBRARY])
    pages = [lambda: pp.box_with_points(by_library, CELL_TYPE, VALUE, metric.label, analysis_id,
                                        f"Cell type {metric.noun} per library. Each point is a library.",
                                        CELL_TYPE_LABEL, color_by=LIBRARY, colors=library_colors, order=order,
                                        legend_title="Library")]
    if CONDITION in cells.columns:
        condition_colors = pp.palette_for(cells[CONDITION].dropna())
        by_condition = _scaled(cell_type_table(cells, [CONDITION]), metric)
        pages.insert(0, lambda: pp.dodged_bars(by_condition, CELL_TYPE, VALUE, CONDITION, condition_colors,
                                               metric.label, analysis_id,
                                               f"Cell type {metric.noun} per condition.", CELL_TYPE_LABEL,
                                               order=order))
        pages += [
            lambda: pp.box_with_points(by_library, CELL_TYPE, VALUE, metric.label, analysis_id,
                                       "Each point is a library.", CELL_TYPE_LABEL, hue=CONDITION,
                                       hue_colors=condition_colors, order=order, legend_title="Condition"),
            lambda: pp.box_with_points(by_library, CONDITION, VALUE, metric.label, analysis_id,
                                       "Each point is a library. Variable y-axis.", "Condition", color_by=LIBRARY,
                                       colors=library_colors, facet=CELL_TYPE, sharey=False,
                                       order=list(condition_colors), legend_title="Library"),
        ]
    return pages


def cell_type_metric_report(metric, filename, cells, donor_cell_type, analysis_id, out_dir):
    """
    One PDF with the pages of a metric: overall per cell type, by condition and library if there are conditions, and
    per donor if donor_cell_type is given.

    :return: the path, or None if there are no cell types
    """
    overall = _scaled(cell_type_table(cells), metric)
    if overall.empty:
        return None
    overall = overall.dropna(subset=[VALUE]).sort_values(VALUE, ascending=False, kind="stable")
    order = list(overall[CELL_TYPE])
    cell_type_colors = pp.palette_for(overall[CELL_TYPE])
    pages = [lambda: pp.sorted_bar_with_median(overall, CELL_TYPE, VALUE, metric.label, cell_type_colors,
                                               analysis_id, f"Overall cell type {metric.noun}.",
                                               xlabel=CELL_TYPE_LABEL, show_median=False)]
    pages += library_pages(cells, metric, analysis_id, order)
    if donor_cell_type is not None:
        donor_colors = pp.donor_palette(donor_cell_type[DONOR])
        pages += donor_cell_type_pages(donor_cell_type, metric, donor_colors, analysis_id)
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.{filename}"))


def cell_type_counts_report(cells, donor_cell_type, analysis_id, out_dir):
    return cell_type_metric_report(COUNTS, "cell_type_counts.pdf", cells, donor_cell_type, analysis_id, out_dir)


def cell_type_median_umis_report(cells, donor_cell_type, analysis_id, out_dir):
    return cell_type_metric_report(MEDIAN_UMIS, "cell_type_median_UMIs.pdf", cells, donor_cell_type, analysis_id,
                                   out_dir)


def cell_type_proportions_report(cells, donor_cell_type, analysis_id, out_dir):
    return cell_type_metric_report(PROPORTIONS, "cell_type_proportions.pdf", cells, donor_cell_type, analysis_id,
                                   out_dir)


def donor_cell_type_distribution_pages(donor_cell_type, analysis_id):
    """
    Stacked cell type percentages per donor, which the donor summary also shows, split into expected and unexpected
    donors.
    """
    donor_colors = pp.donor_palette(donor_cell_type[DONOR])
    cell_type_colors = pp.palette_for(donor_cell_type[CELL_TYPE])
    present = donor_cell_type[donor_cell_type["num_nuclei"] > 0].assign(
        expected=lambda d: d["expected"].map({True: "EXPECTED", False: "UNEXPECTED"}))
    return [lambda: pp.stacked_bars(present, DONOR, CELL_TYPE, "num_nuclei", cell_type_colors, "% of nuclei",
                                    analysis_id, "Donor cell type proportions. Percentages > 3% labeled.", "Donor",
                                    facet="expected")]


def donor_cell_type_distributions_report(cells, donor_cell_type, analysis_id, out_dir):
    """
    Write <analysis_id>.donor_cell_type_distributions.pdf: cell type proportions of each donor, overall, by
    condition and by library.
    """
    pages = donor_cell_type_distribution_pages(donor_cell_type, analysis_id)
    donors = cells.dropna(subset=[DONOR])
    library_colors = pp.palette_for(donors[LIBRARY])
    cell_type_order = list(cell_type_table(cells).sort_values("num_nuclei", ascending=False)[CELL_TYPE])
    if CONDITION in cells.columns:
        condition_colors = pp.palette_for(donors[CONDITION])
        by_condition = cell_type_table(donors, [DONOR, CONDITION])
        by_condition = by_condition[by_condition["num_nuclei"] > 0]
        pages.append(lambda: pp.dodged_bars(by_condition, CELL_TYPE, "fraction_nuclei", CONDITION, condition_colors,
                                            "Fraction of nuclei", analysis_id,
                                            "Donor cell type distributions for each condition.", CELL_TYPE_LABEL,
                                            order=cell_type_order, facet=DONOR))
    by_library = cell_type_table(donors, [DONOR, LIBRARY])
    by_library = by_library[by_library["num_nuclei"] > 0]
    pages.append(lambda: pp.box_with_points(by_library, CELL_TYPE, "fraction_nuclei", "Fraction of nuclei",
                                            analysis_id, "Donor cell type distributions. Each point is a library.",
                                            CELL_TYPE_LABEL, color_by=LIBRARY, colors=library_colors,
                                            facet=DONOR, order=cell_type_order, legend_title="Library"))
    return pp.write_pdf(pages, os.path.join(out_dir, f"{analysis_id}.donor_cell_type_distributions.pdf"))
