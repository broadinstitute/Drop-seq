#!/usr/bin/env python3
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
Summarize QC metrics of Nextflow scRNA-seq experiments, one row per experiment, and optionally write a
standard analysis tear sheet per experiment.

This is a port of buildSummaryTable and makeStandardAnalysisTearSheet from the DropSeq.barnyard.private R package,
adapted to the Nextflow pipeline outputs.  The input artifacts are found by locate_scRNA_artifacts: supply either
DGEs (--dge or --manifest, as for locate_scRNA_artifacts) or its output (--artifacts).

Metrics whose input artifacts are missing are NA.  Relative to the R summary, flowcell, pct_polyATrimmed,
pct_TSOTrimmed, BarcodeIndexPDF and the scPred columns are not produced; scDblFinder_pct_doublets replaces
scPred_doublets.  run_status and run_date come from the most recent Nextflow execution trace.  The tear sheet has no
return to sequencing coverage or sample index rows, and its median UMIs(STAMPS) is the median UMIs of the selected
cells after CBRB rather than the (incorrect) value in the cell selection summary.
"""

import argparse
import io
import os
import posixpath
import re
import sys

import pandas as pd
import yaml

from dropseq.aggregation import logger, add_log_argument, dctLogLevel
import dropseq.aggregation.locate_scRNA_artifacts as locator
from dropseq.aggregation.locate_scRNA_artifacts import NA
from dropseq.util.storage import default_store

SELECTED_CELL_STAT_COLUMNS = ["pct_genic", "pct_exonic", "pct_intronic", "pct_intergenic", "pct_ribosomal", "pct_mt"]

DROPULATION_COLUMNS = {
    # summary stats column: output column
    "pct_all_doublets": "pct_all_doublets",
    "pct_confident_doublets": "pct_confident_doublets",
    "singlets": "singlets",
    "assignable_singlets": "assignable_singlets",
    "reads_per_umi": "reads_per_umi_donor",
    "totalUMIs": "totalUMIs",
    "equitability": "donor_equitability",
    "cell_equitability": "cell_equitability",
}

# output column: artifact manifest field
PATH_COLUMNS = {
    "AlignmentPDF": "alignment_pdf",
    "CBRB_tear_sheet_path": "cbrb_tearsheet",
    "CBRB_html_path": "cbrb_report",
    "SelectionPlot": "cell_selection_pdf",
    "DropulationReport": "dropulation_report",
    "Metacell": "meta_cell",
    "gmg_dge_path": "dge_gmg",
}

SUMMARY_COLUMNS = [
    "library", "exp_date", "reference", "standard_analysis_dir", "run_status", "run_date",
    "total_reads_M", "pct_mapped_reads", "pct_hq_mapped_reads", "selected_cells_total_reads_M",
    "nonselected_cbcs_total_reads_M",
    "pct_genic", "pct_exonic", "pct_intronic", "pct_intergenic", "pct_ribosomal", "pct_mt", "nonselected_cbcs_pct_mt",
    *[f"selected_cells_{c}" for c in SELECTED_CELL_STAT_COLUMNS],
    "n_STAMPs", "median_cell_UMIs", "median_ambient_UMIs", "pct_ambient", "pct_chimeric_UMIs",
    "median_cell_UMIs_cbrb", "cbrb_retained_umi_pct", "reads_per_umi",
    *DROPULATION_COLUMNS.values(),
    "scDblFinder_pct_doublets",
    *PATH_COLUMNS.keys(),
]

# Nextflow task statuses that mean the task's work is done.
TASK_SUCCESS_STATUSES = {"COMPLETED", "CACHED"}
TRACE_DATE_PATTERN = re.compile(r"execution_trace_(\d{4}-\d{2}-\d{2})_")


def _is_na(url):
    return url is None or url == NA


def _na(columns):
    return {column: None for column in columns}


def read_table(url, store, header=True):
    """
    Read the first table in a tab-separated file, skipping '#' lines and blank lines before it.  The table ends at
    the first blank line after it starts, so a Picard metrics histogram that follows is not included.

    :return: DataFrame, with all columns as strings if header is False.
    """
    lines = []
    with store.open_text(url) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if not line.strip():
                if lines:
                    break
                continue
            lines.append(line)
    if not lines:
        raise ValueError(f"No table found in {url}")
    if header:
        return pd.read_csv(io.StringIO("".join(lines)), sep="\t")
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", header=None, dtype=str)


def read_barcodes(url, store):
    """
    :return: set of cell barcodes, one per line, skipping '#' lines.
    """
    with store.open_text(url) as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}


def _require(record, fields, columns, what):
    """
    :return: True if every field of record is present.  Otherwise log a warning and return False.
    """
    missing = [field for field in fields if _is_na(record.get(field))]
    if missing:
        logger.warning(f"{record.get('uei')}: {what} are NA; missing {', '.join(missing)}")
        return False
    return True


def experiment_info(record, store):
    """
    :return: library, exp_date, reference and standard_analysis_dir columns.
    """
    exp_date = None
    if not _is_na(record.get("root_properties")):
        exp_date = (yaml.safe_load(store.read_text(record["root_properties"])) or {}).get("experimentDate")
    return {
        "library": record["uei"],
        "exp_date": exp_date,
        "reference": posixpath.basename(record["alignment_dir"].rstrip("/")),
        "standard_analysis_dir": posixpath.basename(record["standard_analysis_dir"].rstrip("/")),
    }


def run_status(record, store):
    """
    run_status is COMPLETED if every task in the most recent execution trace has a COMPLETED or CACHED attempt,
    otherwise FAILED.  run_date is the start date of that run, from the trace filename.
    """
    columns = ["run_status", "run_date"]
    if not _require(record, ["pipeline_trace"], columns, "run status"):
        return _na(columns)
    url = record["pipeline_trace"]
    match = TRACE_DATE_PATTERN.search(posixpath.basename(url))
    trace = read_table(url, store)
    if trace.empty:
        status = None
    else:
        succeeded = trace["status"].isin(TASK_SUCCESS_STATUSES).groupby(trace["name"]).any()
        status = "COMPLETED" if succeeded.all() else "FAILED"
    return {"run_status": status, "run_date": match.group(1) if match else None}


def read_quality_metrics(record, store):
    columns = ["total_reads_M", "pct_mapped_reads", "pct_hq_mapped_reads"]
    if not _require(record, ["read_quality_metrics"], columns, "read quality metrics"):
        return _na(columns)
    m = read_table(record["read_quality_metrics"], store).iloc[0]
    return {
        "total_reads_M": round(m["totalReads"] / 1e6, 2),
        "pct_mapped_reads": round(m["mappedReads"] / m["totalReads"] * 100, 2),
        "pct_hq_mapped_reads": round(m["hqMappedReads"] / m["totalReads"] * 100, 2),
    }


def pct_intronic_exonic(record, store):
    columns = ["pct_genic", "pct_exonic", "pct_intronic", "pct_intergenic", "pct_ribosomal"]
    if not _require(record, ["frac_intronic_exonic"], columns, "library functional annotation metrics"):
        return _na(columns)
    m = read_table(record["frac_intronic_exonic"], store).iloc[0]
    fractions = {
        "pct_genic": m["PCT_UTR_BASES"] + m["PCT_CODING_BASES"] + m["PCT_INTRONIC_BASES"],
        "pct_exonic": m["PCT_CODING_BASES"] + m["PCT_UTR_BASES"],
        "pct_intronic": m["PCT_INTRONIC_BASES"],
        "pct_intergenic": m["PCT_INTERGENIC_BASES"],
        "pct_ribosomal": m["PCT_RIBOSOMAL_BASES"],
    }
    return {column: round(value * 100, 2) for column, value in fractions.items()}


def selected_cell_functional_stats(record, store, selected=None):
    """
    Functional annotation metrics from the per-cell features of all cell barcodes, split by whether the cell barcode
    was selected.  pct_mt is the read-weighted mean over all cell barcodes, and so the % MT reads in the library.

    :param selected: set of selected cell barcodes, read from selected_cell_barcodes if not given.
    """
    columns = (["pct_mt", "nonselected_cbcs_pct_mt", "nonselected_cbcs_total_reads_M", "selected_cells_total_reads_M"] +
               [f"selected_cells_{c}" for c in SELECTED_CELL_STAT_COLUMNS])
    if not _require(record, ["cell_features", "selected_cell_barcodes"], columns, "selected cell functional stats"):
        return _na(columns)
    if selected is None:
        selected = read_barcodes(record["selected_cell_barcodes"], store)
    features = read_table(record["cell_features"], store)
    features["pct_exonic"] = features["pct_coding"] + features["pct_utr"]
    features["pct_genic"] = features["pct_exonic"] + features["pct_intronic"]
    is_selected = features["cell_barcode"].isin(selected)
    not_selected = features[~is_selected]
    selected_features = features[is_selected]
    result = {
        "pct_mt": (features["pct_mt"] * features["num_reads"]).sum() / features["num_reads"].sum() * 100,
        "nonselected_cbcs_pct_mt": not_selected["pct_mt"].mean() * 100,
        "nonselected_cbcs_total_reads_M": not_selected["num_reads"].sum() / 1e6,
        "selected_cells_total_reads_M": selected_features["num_reads"].sum() / 1e6,
    }
    for column in SELECTED_CELL_STAT_COLUMNS:
        result[f"selected_cells_{column}"] = selected_features[column].mean() * 100
    return {column: result[column] for column in columns}


def median_selected_cbrb_umis(record, store, selected):
    """
    :return: median UMIs of the selected cells after CBRB, from the selected cells DGE summary.
    """
    summary = read_table(record["dge_summary_selected_cells"], store)
    return summary.loc[summary["CELL_BARCODE"].isin(selected), "NUM_TRANSCRIPTS"].median()


def selection_summary_stats(record, store, selected=None):
    """
    Cell selection metrics from the unfiltered DGE summary.  median_cell_UMIs_cbrb is the median UMIs of the selected
    cells after CBRB, from the selected cells DGE summary.

    :param selected: set of selected cell barcodes, read from selected_cell_barcodes if not given.
    """
    columns = ["n_STAMPs", "median_cell_UMIs", "median_ambient_UMIs", "pct_ambient", "reads_per_umi",
               "median_cell_UMIs_cbrb"]
    result = _na(columns)
    if not _require(record, ["selected_cell_barcodes"], columns, "cell selection summary stats"):
        return result
    if selected is None:
        selected = read_barcodes(record["selected_cell_barcodes"], store)
    if _require(record, ["dge_summary_unfiltered"], ["n_STAMPs", "median_cell_UMIs", "median_ambient_UMIs",
                                                     "pct_ambient", "reads_per_umi"], "cell selection summary stats"):
        summary = read_table(record["dge_summary_unfiltered"], store)
        is_stamp = summary["CELL_BARCODE"].isin(selected)
        stamps = summary[is_stamp]
        result["n_STAMPs"] = len(selected)
        result["median_cell_UMIs"] = stamps["NUM_TRANSCRIPTS"].median()
        result["median_ambient_UMIs"] = summary.loc[~is_stamp, "NUM_TRANSCRIPTS"].median()
        result["pct_ambient"] = round(result["median_ambient_UMIs"] / result["median_cell_UMIs"] * 100, 2)
        if _require(record, ["reads_per_cell"], ["reads_per_umi"], "reads per UMI"):
            reads = read_table(record["reads_per_cell"], store, header=False)
            reads = pd.Series(reads[0].astype(int).values, index=reads[1])
            reads_per_umi = stamps["CELL_BARCODE"].map(reads) / stamps["NUM_TRANSCRIPTS"]
            result["reads_per_umi"] = round(reads_per_umi.median(), 2)
    if _require(record, ["dge_summary_selected_cells"], ["median_cell_UMIs_cbrb"], "CBRB median UMIs"):
        result["median_cell_UMIs_cbrb"] = median_selected_cbrb_umis(record, store, selected)
    return result


def pct_chimeric_umis(record, store):
    if not _require(record, ["chimeric_metrics"], ["pct_chimeric_UMIs"], "chimeric metrics"):
        return {"pct_chimeric_UMIs": None}
    m = read_table(record["chimeric_metrics"], store).iloc[0]
    return {"pct_chimeric_UMIs": round(m["NUM_MARKED_UMIS"] / m["NUM_UMIS"] * 100, 3)}


def cbrb_retained_umi_pct(record, store, selected=None):
    """
    % of the UMIs in selected cells that are retained by CBRB, as in the CBRB tear sheet.

    :param selected: set of selected cell barcodes, read from selected_cell_barcodes if not given.
    """
    if not _require(record, ["cbrb_cell_features", "selected_cell_barcodes"], ["cbrb_retained_umi_pct"],
                    "CBRB retained UMIs"):
        return {"cbrb_retained_umi_pct": None}
    if selected is None:
        selected = read_barcodes(record["selected_cell_barcodes"], store)
    features = read_table(record["cbrb_cell_features"], store)
    features = features[features["cell_barcode"].isin(selected)]
    retained = features["num_retained_transcripts"].sum() / features["num_transcripts"].sum()
    return {"cbrb_retained_umi_pct": round(retained * 100, 1)}


def dropulation_summary_stats(record, store):
    columns = list(DROPULATION_COLUMNS.values())
    if _is_na(record.get("dropulation_summary_stats")):
        # absent for experiments without donor assignment, which is not worth a warning
        return _na(columns)
    m = read_table(record["dropulation_summary_stats"], store).iloc[0]
    return {output: m.get(column) for column, output in DROPULATION_COLUMNS.items()}


def scdblfinder_pct_doublets(record, store):
    """
    % of selected cells that scDblFinder calls doublets, from the standard analysis cell metadata.
    """
    if not _require(record, ["selected_cell_metadata"], ["scDblFinder_pct_doublets"], "scDblFinder doublets"):
        return {"scDblFinder_pct_doublets": None}
    metadata = read_table(record["selected_cell_metadata"], store)
    if metadata.empty:
        return {"scDblFinder_pct_doublets": None}
    return {"scDblFinder_pct_doublets": round((metadata["doublet"] == "doublet").mean() * 100, 2)}


def path_columns(record):
    return {column: None if _is_na(record.get(field)) else record[field] for column, field in PATH_COLUMNS.items()}


def summarize_experiment(record, store):
    """
    :param record: dataset dictionary from locate_scRNA_artifacts.
    :return: dictionary of SUMMARY_COLUMNS, in order.  Missing values are None.
    """
    selected = None
    if not _is_na(record.get("selected_cell_barcodes")):
        selected = read_barcodes(record["selected_cell_barcodes"], store)
    result = {}
    result.update(experiment_info(record, store))
    result.update(run_status(record, store))
    result.update(read_quality_metrics(record, store))
    result.update(pct_intronic_exonic(record, store))
    result.update(selected_cell_functional_stats(record, store, selected))
    result.update(selection_summary_stats(record, store, selected))
    result.update(pct_chimeric_umis(record, store))
    result.update(cbrb_retained_umi_pct(record, store, selected))
    result.update(dropulation_summary_stats(record, store))
    result.update(scdblfinder_pct_doublets(record, store))
    result.update(path_columns(record))
    return {column: result[column] for column in SUMMARY_COLUMNS}


def standard_analysis_tearsheet(record, store):
    """
    :param record: dataset dictionary from locate_scRNA_artifacts.
    :return: list of (label, value) rows.
    """
    rows = []
    if _require(record, ["read_quality_metrics"], ["read quality rows"], "tear sheet read quality rows"):
        m = read_table(record["read_quality_metrics"], store).iloc[0]
        rows += [("PF reads", str(m["totalReads"])),
                 ("HQ mapped reads", str(m["hqMappedReads"])),
                 ("% HQ mapped reads", "%.1f%%" % (m["hqMappedReads"] * 100 / m["totalReads"]))]
    chimeric = "NA"
    if not _is_na(record.get("chimeric_metrics")):
        m = read_table(record["chimeric_metrics"], store).iloc[0]
        chimeric = "%.2f%%" % (m["NUM_MARKED_UMIS"] / m["NUM_UMIS"] * 100)
    rows.append(("% chimeric UMIs", chimeric))
    if _require(record, ["cell_selection_report"], ["cell selection rows"], "tear sheet cell selection rows"):
        m = read_table(record["cell_selection_report"], store).iloc[0]
        rows.append(("median Reads/UMI(STAMPS)", _format_value(round(m["medianReadsPerUMI"], 1))))
        if not pd.isna(m["thresholdMessage"]):
            rows.append(("M/O", str(m["thresholdMessage"])))
        rows.append(("STAMPS", _format_value(m["n_STAMPs"])))
        # TODO: UMIs_STAMPs_median in the cell selection summary is wrong for svm_nuclei.
        # Dropseq.cellselection::CallSTAMPsSvmNuclei computes it over every barcode DropSift classifies as a nucleus,
        # before intersecting with the CBRB non-empties, so CBRB-empty barcodes with 0 retained UMIs pull the median
        # down (SMRI_V11_rxn1: 159, vs 194 over the 3752 selected cells, as in the DropSift PDF).  UMIs_STAMPs_min,
        # droplet_cell_occupancy and readsPerUMI in that file have the same problem.  Until that is fixed, compute the
        # median from the selected cells here instead of reporting UMIs_STAMPs_median.
        median_umis = None
        if _require(record, ["dge_summary_selected_cells", "selected_cell_barcodes"], ["median UMIs(STAMPS)"],
                    "tear sheet median UMIs"):
            median_umis = median_selected_cbrb_umis(record, store,
                                                    read_barcodes(record["selected_cell_barcodes"], store))
            median_umis = None if pd.isna(median_umis) else round(median_umis)
        rows.append(("median UMIs(STAMPS)", _format_value(median_umis)))
        if not pd.isna(m["UMIs_ambient"]):
            rows.append(("ambient peak(UMIs)", str(int(m["UMIs_ambient"]))))
    return rows


def _format_value(value):
    """
    Format like R: NA for missing values, and numbers with up to 15 significant digits.
    """
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return NA
    if isinstance(value, float):
        return format(value, ".15g")
    return str(value)


def summarize_datasets(datasets, store=None):
    """
    :return: DataFrame with SUMMARY_COLUMNS and one row per dataset, in input order.
    """
    rows = []
    for record in datasets:
        logger.info(f"Summarizing {record['uei']}")
        rows.append(summarize_experiment(record, store or default_store(record["alignment_dir"])))
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def write_summary(summary, out, remove_empty_columns=False):
    if remove_empty_columns:
        summary = summary.dropna(axis="columns", how="all")
    out.write("\t".join(summary.columns) + "\n")
    for row in summary.itertuples(index=False):
        out.write("\t".join(_format_value(value) for value in row) + "\n")


def tearsheet_filename(uei):
    return f"{uei}.standard_analysis_tearsheet.txt"


def write_tearsheets(datasets, output_dir, store=None):
    """
    Write one tear sheet per dataset to output_dir/<uei>.standard_analysis_tearsheet.txt.
    """
    ueis = [record["uei"] for record in datasets]
    duplicates = sorted({uei for uei in ueis if ueis.count(uei) > 1})
    if duplicates:
        raise ValueError(f"Cannot write one tear sheet per experiment; duplicate uei values: {duplicates}")
    os.makedirs(output_dir, exist_ok=True)
    for record in datasets:
        rows = standard_analysis_tearsheet(record, store or default_store(record["alignment_dir"]))
        with open(os.path.join(output_dir, tearsheet_filename(record["uei"])), "w") as out:
            out.write("label\tvalue\n")
            for label, value in rows:
                out.write(f"{label}\t{value}\n")


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--dge", help="Standard-analysis or village DGE URL.")
    input_group.add_argument("--manifest", type=argparse.FileType('r'),
                             help="YAML manifest with a top-level 'dges' list of 'dge' entries.")
    input_group.add_argument("--artifacts", type=argparse.FileType('r'),
                             help="Artifact manifest written by locate_scRNA_artifacts.")
    parser.add_argument("--summary", "-o", type=argparse.FileType('w'),
                        help="Write the summary, one row per experiment, to this file.  Use - for stdout.")
    parser.add_argument("--tearsheet-dir",
                        help="Write one <uei>.standard_analysis_tearsheet.txt file per experiment to this directory.")
    parser.add_argument("--remove-empty-columns", action="store_true",
                        help="Omit summary columns that are NA for every experiment.")
    options = parser.parse_args(args)
    if options.summary is None and options.tearsheet_dir is None:
        parser.error("at least one of --summary and --tearsheet-dir is required")
    return options


def main(args=None):
    return run(parse_args(args))


def run(options):
    logger.setLevel(dctLogLevel[options.log_level])
    if options.artifacts is not None:
        datasets = locator.load_artifact_manifest(options.artifacts)
        options.artifacts.close()
    else:
        if options.manifest is not None:
            dge_urls = locator.read_dge_manifest(options.manifest)
            options.manifest.close()
        else:
            dge_urls = [options.dge]
        datasets = locator.locate_datasets(dge_urls)["datasets"]
    if options.summary is not None:
        write_summary(summarize_datasets(datasets), options.summary, options.remove_empty_columns)
        if options.summary is not sys.stdout:
            options.summary.close()
    if options.tearsheet_dir is not None:
        write_tearsheets(datasets, options.tearsheet_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
