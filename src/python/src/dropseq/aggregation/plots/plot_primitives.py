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
Page builders for the aggregation QC reports.  Each function draws one 14 x 7 inch page and returns the matplotlib
Figure, which the caller adds to a PdfPages and closes.  Data are long-format DataFrames; colours are dicts from label
to colour.
"""

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

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
MEDIAN_SUBTITLE = "Dashed line is the median (excluding donors with 0 nuclei ascertained)."


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


def page_title(fig, analysis_id, subtitle=None):
    fig.text(0.01, 0.985, analysis_id, ha="left", va="top", fontsize=14, fontweight="bold")
    if subtitle:
        fig.text(0.01, 0.945, subtitle, ha="left", va="top", fontsize=9, color="#444444")


def new_page(analysis_id, subtitle=None, nrows=1, ncols=1, widths=None, height=None, **kwargs):
    """
    :return: Figure and an nrows x ncols array of Axes with a light grid, titled with analysis_id and subtitle
    """
    fig, axes = plt.subplots(nrows, ncols, figsize=(PAGE_SIZE[0], height or PAGE_SIZE[1]), squeeze=False,
                             gridspec_kw={"width_ratios": widths} if widths else None, **kwargs)
    for ax in axes.flat:
        ax.set_axisbelow(True)
        ax.grid(True, color="#e5e5e5", linewidth=0.6)
    page_title(fig, analysis_id, subtitle)
    return fig, axes


def finish_page(fig, legend_handles=None, legend_title=None, rotate_x=True):
    """
    Leave room for the title and put the legend, if any, outside the axes on the right.
    """
    if rotate_x:
        for ax in fig.axes:
            ax.tick_params(axis="x", labelrotation=45)
            for label in ax.get_xticklabels():
                label.set_horizontalalignment("right")
    right = 0.86 if legend_handles else 1.0
    if legend_handles:
        fig.legend(handles=legend_handles, title=legend_title, loc="center left", bbox_to_anchor=(0.87, 0.5),
                   frameon=False, fontsize=8, ncol=1 + len(legend_handles) // 30)
    fig.tight_layout(rect=(0, 0, right, 0.92))
    return fig


def color_handles(colors, labels=None):
    return [Patch(facecolor=colors[label], edgecolor="black", label=label)
            for label in (labels if labels is not None else colors)]


def _median(values, exclude_zero):
    values = pd.Series(values).dropna()
    if exclude_zero:
        values = values[values != 0]
    return values.median() if len(values) else float("nan")


def _color_labels(ax, labels, red):
    for tick in ax.get_xticklabels():
        if tick.get_text() in red:
            tick.set_color(UNEXPECTED_COLOR)


def _bar_labels(ax, bars, values, rotation=0, size=8, red=(), labels=None):
    for i, (bar, value) in enumerate(zip(bars, values)):
        color = UNEXPECTED_COLOR if labels is not None and labels[i] in red else "black"
        ax.annotate(format_value(value), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 2), textcoords="offset points", ha="center", va="bottom", rotation=rotation,
                    fontsize=size, color=color)


def _bars_with_median(ax, df, x, y, colors, red, label_rotation=0, label_size=8, sort=True,
                      median_exclude_zero=True, show_median=True):
    """
    One bar per row of df on ax, sorted by descending value, with value labels and a dashed median line.
    Missing values are drawn as 0 and labelled NA.  Labels of the x values in red are red.
    """
    df = df.assign(_value=df[y].astype(float))
    if sort:
        df = df.sort_values("_value", ascending=False, na_position="last", kind="stable")
    filled = df["_value"].fillna(0)
    bars = ax.bar(range(len(df)), filled, color=[colors[v] for v in df[x]], edgecolor="black", linewidth=0.8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df[x])
    _bar_labels(ax, bars, df["_value"], label_rotation, label_size, red, list(df[x]))
    median = _median(df["_value"], median_exclude_zero) if show_median else float("nan")
    if not math.isnan(median):
        ax.axhline(median, color="black", linestyle="dashed", linewidth=1)
    ax.set_xlim(-0.6, len(df) - 0.4)
    top = filled.max() if len(filled) else 1
    ax.set_ylim(0, (top if top > 0 else 1) * (1.25 if label_rotation else 1.1))
    _color_labels(ax, df[x], red)
    return df


def sorted_bar_with_median(df, x, y, ylabel, colors, analysis_id, subtitle, xlabel="Donor", red=(), legend=True,
                           show_median=True):
    """
    A bar per row, sorted descending, coloured by x with a dashed median line (unless not show_median).  x values in
    red are drawn with red labels and the subtitle says so.
    """
    if red:
        subtitle = f"{subtitle} Red label = unexpected donor."
    fig, axes = new_page(analysis_id, f"{subtitle}\n{MEDIAN_SUBTITLE}".strip() if show_median else subtitle)
    ax = axes[0, 0]
    shown = _bars_with_median(ax, df, x, y, colors, set(red), show_median=show_median)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return finish_page(fig, color_handles(colors, list(shown[x])) if legend else None, xlabel)


def bar_and_condition_box(df, y, ylabel, colors, analysis_id, subtitle, library="library", condition="condition",
                          y_max=None):
    """
    The library page: a bar per library with value labels and a dashed median (left), and a box plot by condition
    with a point per library (right).
    """
    fig, axes = new_page(analysis_id, f"{ylabel} per library. Dashed line is median.\n{subtitle}".strip(), ncols=2,
                         widths=[3, 2])
    ax = axes[0, 0]
    _bars_with_median(ax, df, library, y, colors, set(), sort=False, median_exclude_zero=False)
    ax.set_xlabel("Library")
    ax.set_ylabel(ylabel)
    top = y_max if y_max is not None else max(df[y].max(), 0) * 1.1 or 1
    ax.set_ylim(0, top)
    box_axes = axes[0, 1]
    _draw_box(box_axes, df, condition, y, color_by=library, colors=colors, order=list(dict.fromkeys(df[condition])))
    box_axes.set_xlabel("Condition")
    box_axes.set_ylabel("")
    box_axes.set_ylim(0, top)
    return finish_page(fig, color_handles(colors, list(df[library])), "Library")


def _draw_box(ax, df, x, y, color_by=None, colors=None, hue=None, hue_colors=None, order=None, point_size=7):
    """
    Box plot over x categories with a point per row.  Points are coloured by color_by, or, with hue, boxes are
    dodged by hue and points are black.
    """
    df = df.dropna(subset=[y])
    if df.empty:
        return
    order = order if order is not None else list(dict.fromkeys(df[x]))
    if hue is not None:
        hue_order = list(hue_colors)
        sns.boxplot(data=df, x=x, y=y, hue=hue, order=order, hue_order=hue_order, palette=hue_colors, dodge=True,
                    fliersize=0, linewidth=0.8, ax=ax, legend=False)
        sns.stripplot(data=df, x=x, y=y, hue=hue, order=order, hue_order=hue_order, dodge=True,
                      palette={h: "black" for h in hue_order}, size=point_size - 3, ax=ax, legend=False)
    else:
        sns.boxplot(data=df, x=x, y=y, order=order, color="white", fliersize=0, linewidth=0.8,
                    boxprops={"edgecolor": "gray"}, whiskerprops={"color": "gray"}, capprops={"color": "gray"},
                    medianprops={"color": "gray"}, ax=ax)
        if color_by is not None:
            sns.stripplot(data=df, x=x, y=y, order=order, hue=color_by, palette=colors, hue_order=list(colors),
                          edgecolor="black", linewidth=0.8, size=point_size, jitter=0.15, ax=ax, legend=False)
        else:
            sns.stripplot(data=df, x=x, y=y, order=order, color="black", size=point_size - 3, ax=ax, legend=False)


def box_with_points(df, x, y, ylabel, analysis_id, subtitle, xlabel, color_by=None, colors=None, hue=None,
                    hue_colors=None, facet=None, order=None, sharey=True, legend_title=None, ncols=None):
    """
    Box plot of y by x, with a point per row of df.  Points are coloured by color_by (a column; colors maps its
    values to colours).  With hue, boxes are dodged by hue (hue_colors maps its values to colours) and the points are
    black.  With facet, one panel per value of the facet column, each with its own y axis unless sharey.

    :param order: the x categories in order; default is descending median of y.
    """
    if order is None:
        order = list(df.dropna(subset=[y]).groupby(x)[y].median().sort_values(ascending=False).index)
    panels = [(None, df)] if facet is None else [(f, d) for f, d in df.groupby(facet, sort=True)]
    ncols = ncols or (1 if facet is None else min(len(panels), math.ceil(math.sqrt(len(panels) * 2))))
    nrows = math.ceil(len(panels) / ncols)
    fig, axes = new_page(analysis_id, subtitle, nrows=nrows, ncols=ncols, height=max(PAGE_SIZE[1], 3 * nrows),
                         sharey=sharey and facet is not None)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, (name, panel) in zip(axes.flat, panels):
        panel_order = order if facet is None else [o for o in order if o in set(panel[x])]
        _draw_box(ax, panel, x, y, color_by, colors, hue, hue_colors, panel_order,
                  point_size=7 if facet is None else 5)
        ax.set_xlabel(xlabel if ax in axes[-1] else "")
        ax.set_ylabel(ylabel if ax in axes[:, 0] else "")
        if name is not None:
            ax.set_title(str(name), fontsize=9)
    legend = None
    if hue is not None:
        legend = color_handles(hue_colors)
    elif color_by is not None:
        legend = color_handles(colors, [k for k in colors if k in set(df[color_by])])
    return finish_page(fig, legend, legend_title)


def dodged_bars(df, x, y, hue, hue_colors, ylabel, analysis_id, subtitle, xlabel, order=None, facet=None):
    """
    Bars of y by x, dodged by hue, with value labels.  With facet, one panel per value (sharing x order and y).
    """
    if order is None:
        order = list(df.groupby(x)[y].sum().sort_values(ascending=False).index)
    panels = [(None, df)] if facet is None else list(df.groupby(facet, sort=True))
    ncols = 1 if facet is None else math.ceil(math.sqrt(len(panels) * 2))
    nrows = math.ceil(len(panels) / ncols)
    fig, axes = new_page(analysis_id, subtitle, nrows=nrows, ncols=ncols, height=max(PAGE_SIZE[1], 3 * nrows))
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, (name, panel) in zip(axes.flat, panels):
        sns.barplot(data=panel, x=x, y=y, hue=hue, order=order, hue_order=list(hue_colors), palette=hue_colors,
                    edgecolor="black", linewidth=0.8, dodge=True, ax=ax, legend=False)
        if facet is None:
            for container in ax.containers:
                ax.bar_label(container, labels=[format_value(v) for v in container.datavalues], fontsize=7,
                             padding=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.margins(y=0.12)
        if name is not None:
            ax.set_title(str(name), fontsize=9)
    return finish_page(fig, color_handles(hue_colors), str(hue).replace("_", " ").capitalize())


def stacked_bars(df, x, stack, y, colors, ylabel, analysis_id, subtitle, xlabel, normalize=True, label_min=0.03,
                 facet=None, annotations=None, x_order=None):
    """
    Stacked bars of y by x, stacked by stack.  With normalize each bar is scaled to 100% and the axis is a
    percentage; segments above label_min (a fraction) are labelled with their percentage.  With facet, one panel per
    value, with widths proportional to the number of bars.

    :param annotations: optional dict from x to text drawn above the bar
    """
    df = df.copy()
    panels = [(None, df)] if facet is None else list(df.groupby(facet, sort=True))
    widths = [max(panel[x].nunique(), 1) for _, panel in panels]
    fig, axes = new_page(analysis_id, subtitle, ncols=len(panels), widths=widths)
    stack_order = [s for s in colors if s in set(df[stack])]
    for ax, (name, panel) in zip(axes.flat, panels):
        table = panel.pivot_table(index=x, columns=stack, values=y, aggfunc="sum", fill_value=0)
        table = table.reindex(columns=stack_order, fill_value=0)
        if x_order is not None:
            table = table.reindex([v for v in x_order if v in table.index])
        totals = table.sum(axis=1)
        shares = table.div(totals.replace(0, np.nan), axis=0).fillna(0) if normalize else table
        bottom = np.zeros(len(table))
        for column in stack_order:
            heights = shares[column].to_numpy()
            ax.bar(range(len(table)), heights, bottom=bottom, color=colors[column], edgecolor="black",
                   linewidth=0.8)
            for i, h in enumerate(heights):
                if normalize and h > label_min:
                    ax.text(i, bottom[i] + h / 2, f"{h * 100:.3g}", ha="center", va="center", fontsize=7)
                elif not normalize and h > 0 and h >= label_min * totals.iloc[i]:
                    ax.text(i, bottom[i] + h / 2, format_value(h), ha="center", va="center", fontsize=7)
            bottom += heights
        ax.set_xticks(range(len(table)))
        ax.set_xticklabels(table.index)
        ax.set_xlim(-0.6, len(table) - 0.4)
        ax.set_xlabel(xlabel)
        if annotations:
            for i, key in enumerate(table.index):
                if key in annotations:
                    ax.text(i, bottom[i] * 1.005, annotations[key], ha="center", va="bottom", fontsize=7)
        if normalize:
            ax.set_ylim(0, 1.06 if annotations else 1)
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.set_yticklabels([f"{t * 100:.0f}" for t in np.linspace(0, 1, 6)])
        if name is not None:
            ax.set_title(str(name), fontsize=9)
    axes[0, 0].set_ylabel(ylabel)
    return finish_page(fig, color_handles(colors, stack_order), str(stack).replace("_", " ").capitalize())


def faceted_sorted_bars(df, facet, x, y, colors, ylabel, analysis_id, subtitle, xlabel="Donor", red_by=None,
                        nrows=2):
    """
    One panel per facet value, each with its own y axis: a bar per x value sorted within the panel, rotated value
    labels and the dashed median of the panel (excluding 0 and NA).

    :param red_by: column of df; the x values of rows where it is False have red labels
    """
    panels = list(df.groupby(facet, sort=True))
    nrows = min(nrows, len(panels))
    ncols = math.ceil(len(panels) / nrows)
    unexpected = set(df.loc[~df[red_by].astype(bool), x]) if red_by else set()
    if unexpected:
        subtitle += " Red label = unexpected donor."
    fig, axes = new_page(analysis_id, f"{subtitle}\nDashed line is the median of each panel (excluding 0 and NA). "
                                      "Variable y-axis.", nrows=nrows, ncols=ncols, height=max(PAGE_SIZE[1], 3 * nrows))
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, (name, panel) in zip(axes.flat, panels):
        _bars_with_median(ax, panel, x, y, colors, unexpected, label_rotation=90, label_size=6)
        ax.set_title(str(name), fontsize=9)
        ax.tick_params(axis="x", labelsize=6)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel)
    return finish_page(fig, color_handles(colors), xlabel)


def write_pdf(page_makers, path):
    """
    Write the pages to a PDF, closing each figure after it is saved.

    :param page_makers: callables that each return a Figure, called one at a time to limit memory
    :return: the path, or None if there are no pages (no file is written)
    """
    if not page_makers:
        return None
    with PdfPages(path) as pdf:
        for make_page in page_makers:
            fig = make_page()
            pdf.savefig(fig)
            plt.close(fig)
    logger.info(f"Wrote {path}")
    return path
