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
Page builders for the aggregation QC reports, written with plotnine (a port of R's ggplot2) so that the pages look like
the legacy R reports.  Each builder returns a plotnine plot, which write_pdf draws on one page of a PDF.  Data are
long-format DataFrames.

Colours are part of the look and are the same as in the R reports:
- donors have the 28 colours of PALETTE, assigned to the sorted donors, so a donor has the same colour on every page;
- everything else (cell types, libraries, conditions) has ggplot's default hue palette, assigned to the sorted labels
  independently of the order of the bars.  plotnine's default hue scale gives the same colours as R, which
  test_plot_primitives checks.
"""

import math
import os

# Set before matplotlib, which plotnine draws with, is imported: a headless container has no display.
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt  # only to close the figures after they are written
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from plotnine import (aes, after_scale, element_text, facet_grid, facet_wrap, geom_boxplot, geom_col, geom_hline,
                      geom_point, geom_text, ggplot, labs, position_dodge, position_jitter, position_stack,
                      scale_color_manual, scale_fill_hue, scale_fill_manual, scale_x_discrete, scale_y_continuous,
                      theme, theme_bw)

from dropseq.aggregation import logger

PAGE_SIZE = (14, 7)

# The palette of the legacy donor plots.
PALETTE = [
    "#3478b5", "#aec8e8", "#f37f11", "#f6bc78", "#30a12b", "#99df8a", "#d72729",
    "#f59897", "#9567bd", "#c5b0d5", "#8c564c", "#c49d94", "#e477c2", "#f7b7d3",
    "#7f7f7f", "#c7c7c7", "#d1d11b", "#f0f08d", "#159696", "#8fdbdb", "#d93d6c",
    "#d6819a", "#76ad17", "#7c8f5a", "#1bafe0", "#8ed6ed", "#6729cc", "#9f91eb"
]

UNEXPECTED_COLOR = "red"
# Colour of the value labels of expected and unexpected donors.  The legacy plots had a legend for this; it is dropped.
LABEL_COLORS = {"EXPECTED": "black", "UNEXPECTED": UNEXPECTED_COLOR}
MEDIAN_SUBTITLE = "Dashed line represents the median per donor (excluding donors with 0 nuclei ascertained)."

# Text sizes in points.  R sizes geom_text in mm: 3 mm is 8.5 pt and 2.5 mm is 7 pt.
LABEL_SIZE = 8.5
SMALL_LABEL_SIZE = 7
# Fraction of the width of a two-panel page left of its bars, for the overhang of the first rotated label.
LEFT_MARGIN = 0.06
# Separates the x value from the facet in the x key of the bars of a panel (R's tidytext::reorder_within).
KEY_SEP = "\x1f"


def palette_for(labels):
    """
    :return: dict from each distinct label, in sorted order, to a colour of PALETTE, which is cycled (with a warning)
    if there are more labels than colours
    """
    labels = sorted(set(labels))
    if len(labels) > len(PALETTE):
        logger.warning(f"{len(labels)} labels but only {len(PALETTE)} colours; colours are reused")
    return {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(labels)}


def donor_palette(donors):
    """
    :return: dict from donor to colour, assigned to the sorted donors so a donor has the same colour on every page
    """
    return palette_for(donors)


def format_value(value):
    """
    Value label like R's: whole numbers without decimals, otherwise 4 significant digits (1 decimal from 100 up)
    without trailing zeros, so 0.013 stays 0.013 and 6.7 is not 6.70.
    """
    if pd.isna(value):
        return "NA"
    if float(value).is_integer():
        return f"{int(value):d}"
    return f"{value:.1f}".rstrip("0").rstrip(".") if abs(value) >= 100 else f"{value:.4g}"


def _unique(values):
    return list(dict.fromkeys(values))


def _sorted_category(values):
    """
    :return: values as a categorical whose levels are the sorted distinct values, which is what R's factor does and
    what assigns the colours of the default hue scale, whatever the order of the bars.
    """
    return pd.Categorical(values, categories=sorted(set(values)))


def _ordered_category(values, order):
    return pd.Categorical(values, categories=_unique(order), ordered=True)


def _median(values, exclude_zero):
    values = pd.Series(values).dropna()
    if exclude_zero:
        values = values[values != 0]
    return values.median() if len(values) else float("nan")


def _fill_scale(colors, present, name):
    """
    :param colors: dict from value to colour, or None for the default hue palette
    :param present: the values of the data, so that the legend lists only them
    """
    if colors:
        present = set(present)
        return scale_fill_manual(values=colors, limits=[c for c in colors if c in present], name=name)
    return scale_fill_hue(name=name)


def _theme(width=PAGE_SIZE[0], height=PAGE_SIZE[1], x_size=None):
    """
    theme_bw with the x labels at 45 degrees ending at their tick, as in the R reports.  Do not add a rotation_mode or va:
    plotnine then reserves too little room for the rotated labels and the page cuts them off.
    """
    return theme_bw() + theme(figure_size=(width, height),
                              axis_text_x=element_text(rotation=45, ha="right",
                                                       **({"size": x_size} if x_size else {})))


def _label_gap(top):
    """
    Distance above a bar of its value label, like R's vjust=-0.5.
    """
    return 0.012 * (top if top and top > 0 else 1)


def _bar_page(data, x, y, ylabel, colors, title, subtitle, xlabel, legend_title, legend, show_median, median_exclude_zero,
              red=(), label_size=LABEL_SIZE, y_limits=None):
    """
    The plot of a bar per row of data (in the order of its rows) with its value label and, if show_median, the
    dashed median.  Bars are coloured by x.  The labels of the x values in red are red.
    """
    data = data.reset_index(drop=True)
    value = data[y].astype(float)
    top = y_limits[1] if y_limits else (value.fillna(0).max() if len(value) else 1)
    data = data.assign(_x=_ordered_category(data[x], data[x]), _fill=_sorted_category(data[x]),
                       _height=value.fillna(0), _label=[format_value(v) for v in value],
                       _expected=np.where(data[x].isin(set(red)), "UNEXPECTED", "EXPECTED"))
    data["_ytext"] = data["_height"] + _label_gap(top)
    plot = (ggplot(data, aes("_x", "_height", fill="_fill"))
            + geom_col(color="black")
            + geom_text(aes(y="_ytext", label="_label", color="_expected"), va="bottom", size=label_size)
            + scale_color_manual(values=LABEL_COLORS, guide=None)
            + _fill_scale(colors, data[x], legend_title)
            + scale_y_continuous(limits=y_limits, expand=(0.05, 0, 0.12, 0))
            + labs(title=title, subtitle=subtitle, x=xlabel, y=ylabel)
            + _theme())
    median = _median(value, median_exclude_zero) if show_median else float("nan")
    if not math.isnan(median):
        plot = plot + geom_hline(yintercept=median, linetype="dashed")
    if not legend:
        plot = plot + theme(legend_position="none")
    return plot


def sorted_bar_with_median(df, x, y, ylabel, colors, analysis_id, subtitle, xlabel="Donor", red=(), legend=True,
                           show_median=True):
    """
    A bar per row, sorted descending, coloured by x (with colors, a dict from x to colour, or else the default hue
    palette) with a dashed median line (unless not show_median).  x values in red have red value labels and the
    subtitle says so.
    """
    if red:
        subtitle = f"{subtitle} Red label = unexpected donor."
    data = df.assign(_value=df[y].astype(float)).sort_values("_value", ascending=False, na_position="last",
                                                              kind="stable")
    return _bar_page(data, x, y, ylabel, colors, analysis_id,
                     f"{subtitle}\n{MEDIAN_SUBTITLE}".strip() if show_median else subtitle, xlabel, xlabel, legend,
                     show_median, True, red)


def bar_and_condition_box(df, y, ylabel, colors, analysis_id, subtitle, library="library", condition="condition",
                          y_max=None, width=PAGE_SIZE[0], median=True):
    """
    The library page: a bar per library with value labels and, if median, a dashed median (left), and a box plot by
    condition with a point per library (right).  Libraries are coloured with colors, or the default hue palette.
    """
    top = y_max if y_max is not None else max(df[y].max(), 0) * 1.1 or 1
    left = _bar_page(df, library, y, ylabel, colors, analysis_id, subtitle, "Library", "Library", False, median, False,
                     y_limits=(0, top))
    box = _box_layers(df, condition, y, ylabel="", xlabel="Condition", color_by=library, colors=colors,
                      order=_unique(df[condition]), legend_title="Library", point_size=3.5, y_limits=(0, top))
    # The first rotated library label hangs over the left edge of the page, so the bars get room on the left.
    return (left + theme(plot_margin_left=LEFT_MARGIN) | box) & theme(figure_size=(width, PAGE_SIZE[1]))


def _box_layers(df, x, y, ylabel, xlabel, title=None, subtitle=None, color_by=None, colors=None, hue=None,
                hue_colors=None, box_fill_by=None, box_colors=None, point_color_by=None, facet=None,
                facet_kind="wrap", scales="fixed", nrow=None, order=None, legend_title=None, point_size=3, jitter=0,
                y_limits=None):
    """
    Box plot over x categories with a point per row, one panel per value of the column facet if given.

    The points are filled by color_by.  Or the boxes are dodged by hue and the points are black.  Or the boxes are
    filled by box_fill_by and the points are black, or red where point_color_by is UNEXPECTED.
    """
    data = df.dropna(subset=[y]).copy()
    order = order if order is not None else list(data.groupby(x)[y].mean().sort_values(ascending=False).index)
    data[x] = _ordered_category(data[x], [o for o in order if o in set(data[x])])
    data = data.dropna(subset=[x])
    position = position_jitter(width=jitter, height=0, random_state=1) if jitter else "identity"
    if hue is not None:
        data["_fill"] = _sorted_category(data[hue])
        plot = (ggplot(data, aes(x, y, fill="_fill"))
                + geom_boxplot(outlier_shape=None)
                + geom_point(aes(x=x, y=y, group="_fill"), inherit_aes=False, fill="black", color="black",
                             size=point_size - 1, position=position_dodge(width=0.75))
                + _fill_scale(hue_colors, data[hue], legend_title))
    elif box_fill_by is not None:
        data["_fill"] = _sorted_category(data[box_fill_by])
        data["_expected"] = data[point_color_by] if point_color_by else "EXPECTED"
        plot = (ggplot(data, aes(x, y, fill="_fill"))
                + geom_boxplot(outlier_shape=None)
                + geom_point(aes(x=x, y=y, color="_expected", fill=after_scale("color")), inherit_aes=False,
                             size=point_size - 1, position=position)
                + scale_color_manual(values=LABEL_COLORS, guide=None)
                + _fill_scale(box_colors, data[box_fill_by], legend_title))
    else:
        if color_by is not None:
            data["_fill"] = _sorted_category(data[color_by])
            plot = (ggplot(data, aes(x, y)) + geom_boxplot(color="gray", outlier_shape=None)
                    + geom_point(aes(fill="_fill"), shape="o", color="black", size=point_size, position=position)
                    + _fill_scale(colors, data[color_by], legend_title))
        else:
            plot = (ggplot(data, aes(x, y)) + geom_boxplot(color="gray", outlier_shape=None)
                    + geom_point(color="black", size=point_size - 1, position=position))
    if facet is not None:
        if facet_kind == "grid":
            plot = plot + facet_grid(f". ~ {facet}", scales=scales, space=scales if scales != "fixed" else "fixed")
        else:
            plot = plot + facet_wrap(facet, nrow=nrow, scales=scales)
    if y_limits:
        plot = plot + scale_y_continuous(limits=y_limits, expand=(0.05, 0, 0.1, 0))
    return plot + labs(title=title, subtitle=subtitle, x=xlabel, y=ylabel) + _theme()


def box_with_points(df, x, y, ylabel, analysis_id, subtitle, xlabel, color_by=None, colors=None, hue=None,
                    hue_colors=None, box_fill_by=None, box_colors=None, point_color_by=None, facet=None,
                    facet_kind="wrap", scales="fixed", order=None, legend_title=None, nrow=None, jitter=0):
    """
    Box plot of y by x, with a point per row of df.  Points are coloured by color_by (a column; colors maps its
    values to colours, or is None for the default hue palette).  With hue, boxes are dodged by hue and the points are
    black.  With box_fill_by, boxes are filled by that column (box_colors maps its values to colours) and the points
    are black, or red for the rows where the column point_color_by is UNEXPECTED.  With facet, one panel per value
    of the facet column, sharing or not their axes according to scales (fixed, free, free_x or free_y); facet_kind
    grid puts them in a row with widths proportional to their number of x values.

    :param order: the x categories in order; default is descending mean of y, like R's reorder.
    """
    return _box_layers(df, x, y, ylabel, xlabel, analysis_id, subtitle, color_by, colors, hue, hue_colors,
                       box_fill_by, box_colors, point_color_by, facet, facet_kind, scales, nrow, order, legend_title,
                       point_size=3, jitter=jitter)


def dodged_bars(df, x, y, hue, hue_colors, ylabel, analysis_id, subtitle, xlabel, order=None, facet=None):
    """
    Bars of y by x, dodged by hue (coloured with hue_colors, or the default hue palette if None), with value
    labels.  With facet, one panel per value (sharing x order and axes).
    """
    if order is None:
        order = list(df.groupby(x)[y].sum().sort_values(ascending=False).index)
    data = df.copy()
    data[x] = _ordered_category(data[x], order)
    data = data.dropna(subset=[x])
    data["_fill"] = _sorted_category(data[hue])
    top = data[y].max()
    data["_label"] = [format_value(v) for v in data[y]]
    data["_ytext"] = data[y] + _label_gap(top)
    plot = (ggplot(data, aes(x, y, fill="_fill"))
            + geom_col(position="dodge", color="black")
            + _fill_scale(hue_colors, data[hue], str(hue).replace("_", " ").capitalize())
            + scale_y_continuous(expand=(0.05, 0, 0.12, 0))
            + labs(title=analysis_id, subtitle=subtitle, x=xlabel, y=ylabel)
            + _theme())
    if facet is None:
        plot = plot + geom_text(aes(y="_ytext", label="_label"), position=position_dodge(width=0.9), va="bottom",
                                size=SMALL_LABEL_SIZE)
    else:
        plot = plot + facet_wrap(facet)
    return plot


def stacked_bars(df, x, stack, y, colors, ylabel, analysis_id, subtitle, xlabel, normalize=True, label_min=0.03,
                 facet=None, annotations=None, x_order=None, percent=True, label_segments=True):
    """
    Stacked bars of y by x, stacked by stack (coloured with colors, or the default hue palette if None).  With
    normalize each bar is scaled to a fraction of 1 and, with percent, the axis shows percentages.  Segments above
    label_min (a fraction) are labelled with their percentage, if label_segments.  With facet, one panel per value,
    with widths proportional to the number of bars.

    :param annotations: optional dict from x to text drawn above the bar
    """
    stack_order = [s for s in colors if s in set(df[stack])] if colors else sorted(set(df[stack]))
    panels = [(None, df)] if facet is None else list(df.groupby(facet, sort=True))
    frames = []
    for name, panel in panels:
        table = panel.pivot_table(index=x, columns=stack, values=y, aggfunc="sum", fill_value=0)
        table = table.reindex(columns=stack_order, fill_value=0)
        if x_order is not None:
            table = table.reindex([v for v in x_order if v in table.index])
        totals = table.sum(axis=1)
        shares = table.div(totals.replace(0, np.nan), axis=0).fillna(0) if normalize else table
        long = shares.reset_index().melt(id_vars=x, var_name=stack, value_name="_height")
        long["_total"] = long[x].map(totals)
        long["_panel"] = "" if name is None else str(name)
        frames.append(long)
    data = pd.concat(frames, ignore_index=True)
    data = data[data["_height"] > 0] if normalize else data
    x_levels = _unique(v for f in frames for v in f[x])
    data[x] = _ordered_category(data[x], x_levels)
    data["_fill"] = pd.Categorical(data[stack], categories=stack_order)
    if normalize:
        keep = label_segments & (data["_height"] > label_min)
        data["_label"] = np.where(keep, [f"{h * 100:.3g}" for h in data["_height"]], "")
    else:
        keep = label_segments & (data["_height"] > 0) & (data["_height"] >= label_min * data["_total"])
        data["_label"] = np.where(keep, [format_value(h) for h in data["_height"]], "")
    plot = (ggplot(data, aes(x, "_height", fill="_fill"))
            + geom_col(position=position_stack(), color="black")
            + _fill_scale(colors, data[stack], str(stack).replace("_", " ").capitalize())
            + labs(title=analysis_id, subtitle=subtitle, x=xlabel, y=ylabel)
            + _theme())
    if label_segments:
        plot = plot + geom_text(aes(label="_label"), position=position_stack(vjust=0.5), size=SMALL_LABEL_SIZE + 1)
    upper = 0.1 if annotations else 0.05
    if normalize and percent:
        plot = plot + scale_y_continuous(labels=lambda ticks: [f"{t * 100:.0f}" for t in ticks],
                                         expand=(0.05, 0, upper, 0))
    else:
        plot = plot + scale_y_continuous(expand=(0.05, 0, upper, 0))
    if annotations:
        ann = data.drop_duplicates([x, "_panel"]).copy()
        ann = ann[ann[x].astype(str).isin(set(annotations))]
        ann["_top"] = 1.0 if normalize else ann["_total"]
        ann["_note"] = ann[x].astype(str).map(annotations)
        plot = plot + geom_text(aes(x=x, y="_top", label="_note"), data=ann, inherit_aes=False, va="bottom",
                                size=SMALL_LABEL_SIZE, nudge_y=0.01)
    if facet is not None:
        plot = plot + facet_grid(". ~ _panel", scales="free", space="free")
    return plot


def faceted_sorted_bars(df, facet, x, y, colors, ylabel, analysis_id, subtitle, xlabel="Donor", red_by=None,
                        nrows=2, drop_na=False, median_note=None):
    """
    One panel per facet value, each with its own y axis: a bar per x value sorted within the panel, rotated value
    labels and the dashed median of the panel (excluding 0 and NA).

    :param red_by: column of df; the value labels of the x values of rows where it is False are red
    :param drop_na: if true, rows with no value get no bar, instead of a bar of 0 labelled NA, as in R's median UMIs page
    :param median_note: optional second line of the subtitle.  The note on the red labels, if there are any, is the last.
    """
    if drop_na:
        df = df.dropna(subset=[y])
    unexpected = set(df.loc[~df[red_by].astype(bool), x]) if red_by else set()
    frames, medians = [], []
    for name, panel in df.groupby(facet, sort=True):
        panel = panel.assign(_value=panel[y].astype(float)).sort_values("_value", ascending=False,
                                                                        na_position="last", kind="stable")
        top = panel["_value"].fillna(0).max()
        panel = panel.assign(_height=panel["_value"].fillna(0), _label=[format_value(v) for v in panel["_value"]],
                             _key=[f"{v}{KEY_SEP}{name}" for v in panel[x]],
                             _expected=np.where(panel[x].isin(unexpected), "UNEXPECTED", "EXPECTED"))
        panel["_ytext"] = panel["_height"] + _label_gap(top)
        frames.append(panel)
        medians.append({facet: name, "_median": _median(panel["_value"], True)})
    data = pd.concat(frames, ignore_index=True)
    data["_key"] = _ordered_category(data["_key"], data["_key"])
    data["_fill"] = _sorted_category(data[x])
    medians = pd.DataFrame(medians).dropna()
    panels = len(frames)
    plot = (ggplot(data, aes("_key", "_height", fill="_fill"))
            + geom_col(color="black")
            + geom_text(aes(y="_ytext", label="_label", color="_expected"), va="bottom", angle=90,
                        size=SMALL_LABEL_SIZE)
            + scale_color_manual(values=LABEL_COLORS, guide=None)
            + _fill_scale(colors, data[x], xlabel)
            + scale_x_discrete(labels=lambda keys: [k.split(KEY_SEP)[0] for k in keys])
            + scale_y_continuous(expand=(0.05, 0, 0.2, 0))
            + facet_wrap(facet, nrow=min(nrows, panels), scales="free")
            + labs(title=analysis_id, subtitle="\n".join(
                [subtitle] + ([median_note] if median_note else [])
                + (["Red label = unexpected donor."] if unexpected else [])), x=xlabel, y=ylabel)
            + _theme(x_size=6))
    if len(medians):
        plot = plot + geom_hline(aes(yintercept="_median"), data=medians, linetype="dashed")
    return plot


def write_pdf(page_makers, path):
    """
    Write the pages to a PDF, closing each figure after it is saved.

    :param page_makers: callables that each return a plotnine plot, called one at a time to limit memory
    :return: the path, or None if there are no pages (no file is written)
    """
    if not page_makers:
        return None
    with PdfPages(path) as pdf:
        for make_page in page_makers:
            fig = make_page().draw()
            pdf.savefig(fig)
            plt.close(fig)
    logger.info(f"Wrote {path}")
    return path
