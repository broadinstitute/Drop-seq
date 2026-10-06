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
Locate the artifacts of a Nextflow scRNA-seq workflow run, starting from a DGE in the standard_analysis directory
or in the village directory below it.

The input is either a single DGE URL (--dge) or a YAML manifest (--manifest) of the form:

    dges:
      - dge: gs://bucket/.../experiment1.selected.digital_expression.txt.gz
      - dge: gs://bucket/.../experiment2.donors.digital_expression.txt.gz
        donor: N1

A dgeDefaults dictionary is projected onto each entry that does not set the value itself.  Only dge and donor are
used here, and the other keys of the aggregation manifest (filters, joins, ...) are ignored.  The DGE of each entry is
recorded as user_dge, and its donor as donor (NA if not set).  The command line run fails if two DGEs are from the same
experiment, i.e. have the same uei.

Ancestor stage directories (root, alignment, cbrb, cell_selection) are derived from the DGE path.  Optional
descendants of standard_analysis (mmc, village) are discovered by listing.  The output is a YAML document with a
top-level 'datasets' list containing one record per input DGE, in input order.  Missing artifacts are NA.

DGEs and their summaries are emitted as adjacent pairs named by the DGE type, i.e. the part of the filename between
the uei and .digital_expression.txt.gz: dge_unfiltered (alignment), dge_selected_cells (standard_analysis, type
'selected'), dge_gmg (standard_analysis) and dge_donors (village).  The input DGE is reported in the pair for its
type; any other type (e.g. metagene) gets its own pair after dge_gmg.

cell_metadata is the village cell metadata if there is one, else the standard_analysis cell metadata;
selected_cell_metadata is always the standard_analysis cell metadata, which covers all selected cells.
pipeline_trace is the most recent Nextflow execution trace in <root>/pipeline_info.

MMC outputs are keyed by MMC model (the mmc/<model> directory name), with the files in that directory nested below.
"""

import argparse
import copy
import os
import posixpath
import re
import sys

import yaml

from dropseq.aggregation import logger, add_log_argument, dctLogLevel
from dropseq.util.storage import default_store

NA = "NA"
PROPERTIES_FILE = "properties.yaml"
DGE_SUFFIX = ".digital_expression.txt.gz"
DGE_SUMMARY_SUFFIX = ".digital_expression_summary.txt"
CELL_METADATA_TEMPLATE = "{uei}.cmd.tsv"
REDUCED_GTF_SUFFIX = ".reduced.gtf"
PIPELINE_INFO_DIR = "pipeline_info"
PIPELINE_TRACE_PATTERN = r"^execution_trace_.*\.txt$"

# Relative paths of each stage directory under the Nextflow experiment root.
NF_STAGES = {
    "root": r"^$",
    "alignment": r"^[^/]+$",
    "cbrb": r"^[^/]+/cbrb/[^/]+$",
    "cell_selection": r"^[^/]+/cbrb/[^/]+/cell_selection/[^/]+$",
    "std_analysis": r"^[^/]+/cbrb/[^/]+/cell_selection/[^/]+/standard_analysis$",
    "mmc": r"^[^/]+/cbrb/[^/]+/cell_selection/[^/]+/standard_analysis/mmc/[^/]+$",
    "village": r"^[^/]+/cbrb/[^/]+/cell_selection/[^/]+/standard_analysis/village/[^/]+$",
    "any": r".*",
}


# DGE types whose field name differs from the type in the filename.
DGE_TYPE_FIELD_NAMES = {"selected": "selected_cells"}


def dge_field(dge_type):
    return f"dge_{DGE_TYPE_FIELD_NAMES.get(dge_type, dge_type)}"


def dge_summary_field(dge_type):
    return f"dge_summary_{DGE_TYPE_FIELD_NAMES.get(dge_type, dge_type)}"


# Every dataset record contains these fields, in this order.  A non-standard input DGE type adds its pair after
# dge_summary_gmg.
CANONICAL_FIELDS = [
    # UEI / the DGE and donor given by the user / root
    "uei", "user_dge", "donor", "root_properties", "corrected_barcode_metrics", "pipeline_trace",
    # Alignment
    "alignment_dir", "alignment_properties", "reference", "reduced_gtf", "dge_unfiltered", "dge_summary_unfiltered",
    "frac_intronic_exonic",
    "read_quality_metrics", "reads_per_cell", "cell_features", "chimeric_metrics", "alignment_pdf",
    # CBRB
    "cbrb_dir", "cbrb_properties", "cbrb_cell_features", "cbrb_report", "cbrb_tearsheet",
    "cbrb_parameter_estimation_pdf",
    # Cell selection
    "cell_selection_properties", "selected_cell_barcodes", "cell_selection_report", "cell_selection_pdf",
    # Standard analysis
    "standard_analysis_dir", "standard_analysis_properties", "dge_selected_cells", "dge_summary_selected_cells",
    "dge_gmg", "dge_summary_gmg", "cell_metadata", "selected_cell_metadata",
    # MMC: {model: {MMC_FIELDS}}, or NA
    "mmc",
    # Village / donor assignment
    "village_properties", "dge_donors", "dge_summary_donors", "dropulation_summary_stats", "donor_cell_map",
    "donor_assignments", "meta_cell", "dropulation_report", "dropulation_tearsheet",
]

# (field, stage, filename template) for scalar artifacts that live directly in a stage directory.
STAGE_FILES = [
    ("root_properties", "root", PROPERTIES_FILE),
    ("corrected_barcode_metrics", "root", "{uei}.corrected_barcode_metrics"),
    ("alignment_properties", "alignment", PROPERTIES_FILE),
    ("dge_unfiltered", "alignment", "{uei}" + DGE_SUFFIX),
    ("dge_summary_unfiltered", "alignment", "{uei}" + DGE_SUMMARY_SUFFIX),
    ("frac_intronic_exonic", "alignment", "{uei}.fracIntronicExonic.txt"),
    ("read_quality_metrics", "alignment", "{uei}.ReadQualityMetrics.txt"),
    ("reads_per_cell", "alignment", "{uei}.numReads_perCell.txt.gz"),
    ("cell_features", "alignment", "{uei}.cell_features.txt"),
    ("chimeric_metrics", "alignment", "{uei}.chimeric_read_metrics"),
    ("alignment_pdf", "alignment", "{uei}.alignment_summary.pdf"),
    ("cbrb_properties", "cbrb", PROPERTIES_FILE),
    ("cbrb_cell_features", "cbrb", "{uei}.cbrb.cell_features.txt"),
    ("cbrb_report", "cbrb", "{uei}_report.html"),
    ("cbrb_tearsheet", "cbrb", "{uei}.cbrb_tearsheet.pdf"),
    ("cbrb_parameter_estimation_pdf", "cbrb", "{uei}.svm_cbrb_parameter_estimation.pdf"),
    ("cell_selection_properties", "cell_selection", PROPERTIES_FILE),
    ("selected_cell_barcodes", "cell_selection", "{uei}.selectedCellBarcodes.txt"),
    ("cell_selection_report", "cell_selection", "{uei}.cell_selection_assignments_summary.txt"),
    ("cell_selection_pdf", "cell_selection", "{uei}.cell_selection_assignments.pdf"),
    ("standard_analysis_properties", "std_analysis", PROPERTIES_FILE),
    ("dge_selected_cells", "std_analysis", "{uei}.selected" + DGE_SUFFIX),
    ("dge_summary_selected_cells", "std_analysis", "{uei}.selected" + DGE_SUMMARY_SUFFIX),
    ("dge_gmg", "std_analysis", "{uei}.gmg" + DGE_SUFFIX),
    ("dge_summary_gmg", "std_analysis", "{uei}.gmg" + DGE_SUMMARY_SUFFIX),
    ("selected_cell_metadata", "std_analysis", CELL_METADATA_TEMPLATE),
    ("village_properties", "village", PROPERTIES_FILE),
    ("dge_donors", "village", "{uei}.donors" + DGE_SUFFIX),
    ("dge_summary_donors", "village", "{uei}.donors" + DGE_SUMMARY_SUFFIX),
    ("dropulation_summary_stats", "village", "{uei}.dropulation_summary_stats.txt"),
    ("donor_cell_map", "village", "{uei}.donor_cell_map.txt"),
    ("donor_assignments", "village", "{uei}.donor_assignments.txt"),
    ("meta_cell", "village", "{uei}.donors.metacells.txt"),
    ("dropulation_report", "village", "{uei}.dropulation_report.pdf"),
    ("dropulation_tearsheet", "village", "{uei}.dropulation_tearsheet.pdf"),
]

# (field, filename template) for the artifacts in each mmc/<model> directory.
MMC_FILES = [
    ("mmc_properties", PROPERTIES_FILE),
    ("mmc_annotations", "{uei}.csv"),
    ("mmc_json", "{uei}.json"),
    ("mmc_cell_type_counts", "{uei}.cell_type_counts.tsv"),
]


def _parent(url):
    return posixpath.dirname(url.rstrip("/"))


def _basename(url):
    return posixpath.basename(url.rstrip("/"))


def derive_stage_dirs(dge_url):
    """
    Derive the ancestor stage directories of a DGE from its path.  The DGE may be in the standard_analysis
    directory or in a village directory below it.  No I/O is performed.

    :return: dict of stage name to directory URL (no trailing slash) for root, alignment, cbrb, cell_selection
    and std_analysis, plus village if the DGE is in a village directory.
    """
    dge_dir = _parent(dge_url)
    village = None
    if _basename(_parent(dge_dir)) == "village":
        village = dge_dir
        std_analysis = _parent(_parent(dge_dir))
    else:
        std_analysis = dge_dir
    cell_selection = _parent(std_analysis)
    cbrb = _parent(_parent(cell_selection))
    alignment = _parent(_parent(cbrb))
    root = _parent(alignment)
    if (_basename(std_analysis) != "standard_analysis" or
            _basename(_parent(cell_selection)) != "cell_selection" or
            _basename(_parent(cbrb)) != "cbrb"):
        raise ValueError(f"DGE is not in a Nextflow <reference>/cbrb/<x>/cell_selection/<y>/standard_analysis "
                         f"or standard_analysis/village/<z> directory: {dge_url}")
    stage_dirs = {
        "root": root,
        "alignment": alignment,
        "cbrb": cbrb,
        "cell_selection": cell_selection,
        "std_analysis": std_analysis,
    }
    if village is not None:
        stage_dirs["village"] = village
    for stage, stage_dir in stage_dirs.items():
        relative = stage_dir[len(root):].lstrip("/")
        if not re.match(NF_STAGES[stage], relative):
            raise ValueError(f"Derived {stage} directory {stage_dir} does not match the Nextflow layout "
                             f"for DGE {dge_url}")
    return stage_dirs


def locate_reference(alignment_properties, store=None):
    """
    Find the reference of an alignment, and the reduced GTF in the reference directory.

    :param alignment_properties: the alignment properties.yaml, or NA.
    :return: (reference, reduced_gtf), where reference is the 'reference' of the alignment properties, and
    reduced_gtf is the one *.reduced.gtf file in its directory.  Each is NA if it cannot be found.
    :raise ValueError: if there is more than one reduced GTF.
    """
    if alignment_properties == NA:
        return NA, NA
    if store is None:
        store = default_store(alignment_properties)
    reference = (yaml.safe_load(store.read_text(alignment_properties)) or {}).get("reference")
    if not reference:
        logger.warning(f"No reference in {alignment_properties}")
        return NA, NA
    reference = str(reference)
    reference_dir = _parent(reference)
    gtfs = sorted(f for f in default_store(reference_dir).list_dir(reference_dir)[0] if f.endswith(REDUCED_GTF_SUFFIX))
    if len(gtfs) > 1:
        raise ValueError(f"Expected one {REDUCED_GTF_SUFFIX} file in {reference_dir}, found {gtfs}")
    if not gtfs:
        logger.warning(f"No {REDUCED_GTF_SUFFIX} file in {reference_dir}")
        return reference, NA
    return reference, posixpath.join(reference_dir, gtfs[0])


def _resolve(stage_dir, listing, filename):
    return posixpath.join(stage_dir, filename) if filename in listing else NA


def locate_scRNA_artifacts(dge_url, store=None, donor=None):
    """
    Resolve the artifacts related to one DGE in a standard_analysis or village directory.

    :param dge_url: the DGE, named <uei>.<type>.digital_expression.txt.gz.  It is returned unchanged as
    dge_<type>, and as user_dge.
    :param store: GcsStore or LocalStore.  Chosen from the URL scheme if not given.
    :param donor: the donor of a library that is not a village, returned as donor.  NA if not given.
    :return: dataset dictionary with every field in CANONICAL_FIELDS, plus dge_<type> and dge_summary_<type> for a
    non-standard type (inserted after dge_summary_gmg).  Missing artifacts are NA.
    """
    if store is None:
        store = default_store(dge_url)
    dge_name = _basename(dge_url)
    if not dge_name.endswith(DGE_SUFFIX):
        raise ValueError(f"DGE filename does not end with {DGE_SUFFIX}: {dge_url}")
    stage_dirs = derive_stage_dirs(dge_url)
    dge_stage = "village" if "village" in stage_dirs else "std_analysis"

    listings = {stage: store.list_dir(stage_dir) for stage, stage_dir in stage_dirs.items()}
    std_files, std_subdirs = listings["std_analysis"]
    if dge_name not in listings[dge_stage][0]:
        raise FileNotFoundError(f"DGE does not exist: {dge_url}")
    if PROPERTIES_FILE not in listings["root"][0]:
        raise FileNotFoundError(f"Root {PROPERTIES_FILE} not found for DGE {dge_url}")
    root_properties = yaml.safe_load(store.read_text(posixpath.join(stage_dirs["root"], PROPERTIES_FILE))) or {}
    uei = root_properties.get("library")
    if not uei:
        raise ValueError(f"'library' not found in root {PROPERTIES_FILE} for DGE {dge_url}")

    # village is optional, at most one
    village_subdirs = []
    if "village" in std_subdirs:
        village_subdirs = store.list_dir(posixpath.join(stage_dirs["std_analysis"], "village"))[1]
    if len(village_subdirs) > 1:
        raise ValueError(f"Expected at most one village analysis, found {village_subdirs} for DGE {dge_url}")
    if village_subdirs and "village" not in stage_dirs:
        stage_dirs["village"] = posixpath.join(stage_dirs["std_analysis"], "village", village_subdirs[0])
        listings["village"] = store.list_dir(stage_dirs["village"])

    # the DGE type is the part of the filename between the uei and the DGE suffix, e.g. 'selected' or 'donors'
    dge_prefix = dge_name[:-len(DGE_SUFFIX)]
    if not dge_prefix.startswith(uei + ".") or len(dge_prefix) == len(uei) + 1:
        raise ValueError(f"DGE filename is not of the form {uei}.<type>{DGE_SUFFIX}: {dge_url}")
    dge_type = re.sub(r"\W", "_", dge_prefix[len(uei) + 1:])

    fields = list(CANONICAL_FIELDS)
    if dge_field(dge_type) not in fields:
        insert_at = fields.index(dge_summary_field("gmg")) + 1
        fields[insert_at:insert_at] = [dge_field(dge_type), dge_summary_field(dge_type)]
    record = {field: NA for field in fields}
    record["uei"] = uei
    record["user_dge"] = dge_url
    record["donor"] = NA if donor is None else str(donor)
    for field, stage, template in STAGE_FILES:
        if stage in stage_dirs:
            record[field] = _resolve(stage_dirs[stage], listings[stage][0], template.format(uei=uei))
    record["reference"], record["reduced_gtf"] = locate_reference(record["alignment_properties"], store)
    record["alignment_dir"] = stage_dirs["alignment"]
    record["cbrb_dir"] = stage_dirs["cbrb"]
    record["standard_analysis_dir"] = stage_dirs["std_analysis"]

    # trace files are named with the run start time, so the greatest name is the most recent run
    if PIPELINE_INFO_DIR in listings["root"][1]:
        pipeline_info = posixpath.join(stage_dirs["root"], PIPELINE_INFO_DIR)
        traces = sorted(f for f in store.list_dir(pipeline_info)[0] if re.match(PIPELINE_TRACE_PATTERN, f))
        if traces:
            record["pipeline_trace"] = posixpath.join(pipeline_info, traces[-1])

    # the input DGE is reported exactly as given, paired with its summary
    record[dge_field(dge_type)] = dge_url
    record[dge_summary_field(dge_type)] = _resolve(stage_dirs[dge_stage], listings[dge_stage][0],
                                                   dge_prefix + DGE_SUMMARY_SUFFIX)

    # prefer the village cell metadata, which includes donor information
    cell_metadata_name = CELL_METADATA_TEMPLATE.format(uei=uei)
    record["cell_metadata"] = NA
    if "village" in stage_dirs:
        record["cell_metadata"] = _resolve(stage_dirs["village"], listings["village"][0], cell_metadata_name)
    if record["cell_metadata"] == NA:
        record["cell_metadata"] = _resolve(stage_dirs["std_analysis"], std_files, cell_metadata_name)

    # mmc is optional, possibly many, keyed by model
    if "mmc" in std_subdirs:
        mmc_root = posixpath.join(stage_dirs["std_analysis"], "mmc")
        mmc = {}
        for model in store.list_dir(mmc_root)[1]:
            mmc_dir = posixpath.join(mmc_root, model)
            mmc_files = store.list_dir(mmc_dir)[0]
            mmc[model] = {field: _resolve(mmc_dir, mmc_files, template.format(uei=uei))
                          for field, template in MMC_FILES}
        if mmc:
            record["mmc"] = mmc
    return record


def read_dge_entries(file):
    """
    :param file: open file containing a YAML document with a top-level 'dges' list of {dge: URL} entries, and
    optionally a dgeDefaults dictionary, which is projected onto each entry that does not set the value itself.
    :return: list of entry dictionaries, in manifest order.
    """
    manifest = yaml.safe_load(file)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("dges"), list):
        raise ValueError("DGE manifest must have a top-level 'dges' list")
    defaults = manifest.get("dgeDefaults") or {}
    if not isinstance(defaults, dict) or "dge" in defaults:
        raise ValueError("dgeDefaults must be a dictionary without 'dge'")
    entries = []
    for entry in manifest["dges"]:
        if not isinstance(entry, dict) or "dge" not in entry:
            raise ValueError(f"DGE manifest entry has no 'dge' key: {entry}")
        entries.append({**copy.deepcopy(defaults), **entry})
    return entries


def read_dge_manifest(file):
    """
    :param file: open file containing a YAML document with a top-level 'dges' list of {dge: URL} entries.
    :return: list of DGE URLs, in manifest order.
    """
    return [entry["dge"] for entry in read_dge_entries(file)]


def locate_datasets(dge_urls, store=None, donors=None):
    """
    :param donors: donor of each DGE, or None for no donors.  An element may be None for no donor.
    :return: {'datasets': [...]} with one dataset dictionary per DGE, in input order.
    """
    if donors is not None and len(donors) != len(dge_urls):
        raise ValueError(f"{len(donors)} donors for {len(dge_urls)} DGEs")
    datasets = []
    for i, dge_url in enumerate(dge_urls):
        logger.info(f"Locating artifacts for {dge_url}")
        dge_store = store if store is not None else default_store(dge_url)
        datasets.append(locate_scRNA_artifacts(dge_url, dge_store, None if donors is None else donors[i]))
    return {"datasets": datasets}


def check_unique_ueis(result):
    """
    :raise ValueError: if two datasets have the same uei, i.e. two DGEs are from the same experiment.
    """
    ueis = [dataset["uei"] for dataset in result["datasets"]]
    duplicates = sorted({uei for uei in ueis if ueis.count(uei) > 1})
    if duplicates:
        raise ValueError(f"More than one DGE for the same experiment; duplicate uei values: {duplicates}")


def _dump_yaml(result, out):
    yaml.safe_dump(result, out, sort_keys=False, default_flow_style=False, width=sys.maxsize)


def write_artifact_manifest(result, out):
    """
    Write all dataset records to one YAML document.
    """
    _dump_yaml(result, out)


def write_artifact_manifest_dir(result, output_dir):
    """
    Write one YAML document per dataset to output_dir/<uei>.yaml.  Each has a 'datasets' list of length 1.
    """
    ueis = [dataset["uei"] for dataset in result["datasets"]]
    duplicates = sorted({uei for uei in ueis if ueis.count(uei) > 1})
    if duplicates:
        raise ValueError(f"Cannot write one file per experiment; duplicate uei values: {duplicates}")
    os.makedirs(output_dir, exist_ok=True)
    for dataset in result["datasets"]:
        with open(os.path.join(output_dir, dataset["uei"] + ".yaml"), "w") as out:
            _dump_yaml({"datasets": [dataset]}, out)


def load_artifact_manifest(file):
    """
    Parse a single- or multi-dataset artifact manifest.

    :param file: path or open file.
    :return: list of dataset dictionaries.
    """
    if isinstance(file, (str, os.PathLike)):
        with open(file) as f:
            manifest = yaml.safe_load(f)
    else:
        manifest = yaml.safe_load(file)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("datasets"), list):
        raise ValueError("Artifact manifest must have a top-level 'datasets' list")
    return manifest["datasets"]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--dge", help="Standard-analysis DGE URL.")
    parser.add_argument("--donor", help="Donor of the --dge DGE, for a library that is not a village.")
    input_group.add_argument("--manifest", type=argparse.FileType('r'),
                             help="YAML manifest with a top-level 'dges' list of 'dge' entries.")
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument("--output", "-o", type=argparse.FileType('w'), default=sys.stdout,
                              help="Write all datasets to this YAML file.  Default: stdout")
    output_group.add_argument("--output-dir",
                              help="Write one <uei>.yaml file per dataset to this directory.")
    options = parser.parse_args(args)
    if options.donor is not None and options.dge is None:
        parser.error("--donor can only be used with --dge; use a donor entry in the manifest")
    return options


def main(args=None):
    return run(parse_args(args))


def run(options):
    logger.setLevel(dctLogLevel[options.log_level])
    if options.manifest is not None:
        entries = read_dge_entries(options.manifest)
        options.manifest.close()
    else:
        entries = [{"dge": options.dge, "donor": options.donor}]
    result = locate_datasets([entry["dge"] for entry in entries], donors=[entry.get("donor") for entry in entries])
    check_unique_ueis(result)
    if options.output_dir is not None:
        write_artifact_manifest_dir(result, options.output_dir)
    else:
        write_artifact_manifest(result, options.output)
        if options.output is not sys.stdout:
            options.output.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
