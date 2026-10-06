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
Resolve the filters and joins of the manifest of the scRNA aggregation workflow, which has the same schema as the
manifest of Zamboni's LaunchScRnaAggregation:

    dgeDefaults:             # values projected onto each dge that does not set them
      filters: ...
    dges:
      - dge: gs://.../<uei>.donors.digital_expression.txt.gz   # required
        donor: N1            # handled by locate_scRNA_artifacts, which records it as donor in the artifacts
        filters:             # keyed by cell metadata column
          doublet: {exclude: doublet}
          pct_mt: {max: 0.1}      # include, exclude, includeFile, excludeFile, min and max are allowed
        joins:               # secondary files joined on to the cell metadata, in order
          - {joinFile: extra.tsv, leftColumn: cell_barcode, joinColumn: cell}

The library ID is the uei of the DGE, so libraryId is not supported.  The artifacts of each DGE are found by matching
it with user_dge in the per-library artifact files written by locate_scRNA_artifacts --output-dir.  For each DGE,
this writes the arguments for join_and_filter_tsv that apply the donor (--set donor), joins and filters, one argument
per line, and a header-ed tab-separated samplesheet, libraries.tsv, with one row per library:

    library_id  artifact_file  dge  cell_metadata  mmc_model  mmc_annotations  reduced_gtf  args_file

dge is the donor DGE if the library has one, otherwise the selected cells DGE.  mmc_model is --mmc-model, or the first
model of the library, and mmc_annotations is the MapMyCells csv of that model.  Missing values are NA.

The library IDs must differ from the analysis identifier, all libraries must have the same reference (the reference of
the artifacts, or the name of the alignment directory if that is NA) unless --reference is given, and all libraries
must have the same reduced GTF.  Every library must have a DGE, cell metadata and a reduced GTF.

This is a stand-in for the filter and join step of the Nextflow workflow, and goes away with it.
"""

import argparse
import collections
import copy
import os
import posixpath
import sys

import yaml

from dropseq.aggregation import logger, add_log_argument
from dropseq.aggregation.locate_scRNA_artifacts import NA, load_artifact_manifest

DONOR = "donor"
DGE_KEYS = {"dge", "donor", "filters", "joins"}
FILTER_KEYS = {"includeFile", "excludeFile", "include", "exclude", "min", "max"}
JOIN_KEYS = {"joinFile", "leftColumn", "joinColumn"}
LIBRARIES_FILE = "libraries.tsv"
LIBRARIES_COLUMNS = ["library_id", "artifact_file", "dge", "cell_metadata", "mmc_model", "mmc_annotations",
                     "reduced_gtf", "args_file"]
ARGS_FILE_TEMPLATE = "{library_id}.join_and_filter_args.txt"


Library = collections.namedtuple("Library", ["library_id", "artifact_file", "args", "dge", "cell_metadata", "mmc_model",
                                             "mmc_annotations", "reduced_gtf"])


def force_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def load_manifest(file):
    """
    :param file: path or open file.
    """
    if isinstance(file, (str, os.PathLike)):
        with open(file) as f:
            manifest = yaml.safe_load(f)
    else:
        manifest = yaml.safe_load(file)
    if not isinstance(manifest, dict) or "dges" not in manifest:
        raise ValueError("Manifest must have a top-level 'dges' key")
    return manifest


def _check_keys(what, mapping, allowed):
    if not isinstance(mapping, dict):
        raise ValueError(f"{what} must be a dictionary, not {mapping!r}")
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        raise ValueError(f"Unknown key(s) in {what}: {unknown}.  Allowed: {sorted(allowed)}")


def resolve_dges(manifest):
    """
    Project dgeDefaults onto each dge that does not set a value explicitly.

    :return: list of dge dictionaries
    """
    defaults = manifest.get("dgeDefaults") or {}
    _check_keys("dgeDefaults", defaults, DGE_KEYS - {"dge"})
    dges = force_list(manifest["dges"])
    if not dges:
        raise ValueError("At least one DGE must be specified in the manifest")
    resolved = []
    for entry in dges:
        _check_keys("dge entry", entry, DGE_KEYS)
        if "dge" not in entry:
            raise ValueError(f"dge entry has no 'dge' key: {entry}")
        resolved.append({**copy.deepcopy(defaults), **entry})
    return resolved


def _filters_by_column(filters):
    """
    :return: dict of column to filter spec.  A list of dictionaries, each keyed by column, is also accepted.
    """
    if filters is None:
        return {}
    if isinstance(filters, dict):
        return filters
    by_column = {}
    for item in force_list(filters):
        if not isinstance(item, dict):
            raise ValueError(f"filters must be a dictionary keyed by column, not {item!r}")
        for column, spec in item.items():
            if column in by_column:
                raise ValueError(f"More than one filter for column {column}")
            by_column[column] = spec
    return by_column


def join_and_filter_args(dge):
    """
    :param dge: dge dictionary from resolve_dges
    :return: list of join_and_filter_tsv arguments that apply the joins and filters.
    """
    args = []
    for join in force_list(dge.get("joins")):
        _check_keys("join", join, JOIN_KEYS)
        missing = sorted(JOIN_KEYS - set(join))
        if missing:
            raise ValueError(f"join {join} is missing {missing}")
        args += ["--join", str(join["joinFile"]), str(join["leftColumn"]), str(join["joinColumn"])]
    for column, spec in _filters_by_column(dge.get("filters")).items():
        _check_keys(f"filter on {column}", spec, FILTER_KEYS)
        if "includeFile" in spec:
            args += ["--include-file", column, str(spec["includeFile"])]
        if "excludeFile" in spec:
            args += ["--exclude-file", column, str(spec["excludeFile"])]
        if "include" in spec:
            args += ["--include", column, *[str(v) for v in force_list(spec["include"])]]
        if "exclude" in spec:
            args += ["--exclude", column, *[str(v) for v in force_list(spec["exclude"])]]
        if "min" in spec:
            args += ["--min", column, str(spec["min"])]
        if "max" in spec:
            args += ["--max", column, str(spec["max"])]
    return args


def load_artifacts_by_dge(artifacts_dir):
    """
    :return: dict of user_dge to (artifact file, dataset), for every artifact file in the directory.
    """
    by_dge = {}
    for name in sorted(os.listdir(artifacts_dir)):
        if not name.endswith(".yaml"):
            continue
        path = os.path.join(artifacts_dir, name)
        datasets = load_artifact_manifest(path)
        if len(datasets) != 1:
            raise ValueError(f"Expected one dataset in {path}, found {len(datasets)}")
        if datasets[0].get("user_dge", NA) == NA:
            raise ValueError(f"{path} has no user_dge.  Was it written by locate_scRNA_artifacts?")
        by_dge[datasets[0]["user_dge"]] = (path, datasets[0])
    return by_dge


def _mmc(dataset, mmc_model):
    """
    :return: tuple of MapMyCells model and annotations csv.  The model is mmc_model, or the first model of the dataset
    if that is None.  The annotations are NA if the dataset has no such model.
    """
    mmc = dataset.get("mmc")
    mmc = mmc if isinstance(mmc, dict) else {}
    model = mmc_model or next(iter(mmc), NA)
    return model, (mmc.get(model) or {}).get("mmc_annotations", NA)


def resolve(manifest, artifacts_dir, analysis_id, reference=None, mmc_model=None):
    """
    :param manifest: parsed manifest
    :param artifacts_dir: directory of per-library artifact files written by locate_scRNA_artifacts --output-dir
    :param analysis_id: the analysis identifier, which no library ID (uei) may equal
    :param reference: the reference of every library.  If not given, the libraries must all have the same one.
    :param mmc_model: MapMyCells model.  If not given, the first model of each library.
    :return: list of Library, in manifest order.  The args are those for join_and_filter_tsv: a donor from the
    manifest overrides the donor column of the cell metadata.
    """
    by_dge = load_artifacts_by_dge(artifacts_dir)
    resolved = []
    references = {}
    for dge in resolve_dges(manifest):
        if dge["dge"] not in by_dge:
            raise ValueError(f"No artifacts found in {artifacts_dir} for DGE {dge['dge']}")
        artifact_file, dataset = by_dge[dge["dge"]]
        library_id = dataset["uei"]
        library_reference = dataset.get("reference", NA)
        if library_reference == NA:
            library_reference = posixpath.basename(dataset["alignment_dir"].rstrip("/"))
        references.setdefault(library_reference, []).append(library_id)
        args = join_and_filter_args(dge)
        donor = dataset.get(DONOR, NA)
        if donor != NA:
            args = ["--set", "donor", str(donor)] + args
        dge_file = dataset.get("dge_donors", NA)
        if dge_file == NA:
            dge_file = dataset.get("dge_selected_cells", NA)
        cell_metadata = dataset.get("cell_metadata", NA)
        if NA in (dge_file, cell_metadata):
            raise ValueError(f"No DGE or cell metadata found for {library_id}")
        model, annotations = _mmc(dataset, mmc_model)
        resolved.append(Library(library_id, artifact_file, args, dge_file, cell_metadata, model, annotations,
                                dataset.get("reduced_gtf", NA)))
    library_ids = [library.library_id for library in resolved]
    if len(set(library_ids)) != len(library_ids):
        raise ValueError(f"Library IDs must be unique: {library_ids}")
    if analysis_id in library_ids:
        raise ValueError(f"A library ID cannot be the same as the analysis identifier: {analysis_id}")
    if reference is None and len(references) > 1:
        raise ValueError(f"All DGEs must have the same reference.  Found {references}.  "
                         f"Use --reference to override.")
    reduced_gtfs = {library.reduced_gtf for library in resolved}
    if NA in reduced_gtfs:
        raise ValueError(f"No reduced GTF found for {[l.library_id for l in resolved if l.reduced_gtf == NA]}")
    if len(reduced_gtfs) > 1:
        raise ValueError(f"All libraries must have the same reduced GTF.  Found {sorted(reduced_gtfs)}")
    return resolved


def write_resolved(resolved, output_dir):
    """
    Write output_dir/libraries.tsv, and one arguments file per library with one join_and_filter_tsv argument per line.
    """
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, LIBRARIES_FILE), "w") as libraries:
        libraries.write("\t".join(LIBRARIES_COLUMNS) + "\n")
        for library in resolved:
            args_file = os.path.join(output_dir, ARGS_FILE_TEMPLATE.format(library_id=library.library_id))
            with open(args_file, "w") as out:
                out.writelines(arg + "\n" for arg in library.args)
            row = [library.library_id, library.artifact_file, library.dge, library.cell_metadata, library.mmc_model,
                   library.mmc_annotations, library.reduced_gtf, args_file]
            libraries.write("\t".join(row) + "\n")


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_log_argument(parser)
    parser.add_argument("--manifest", "-m", required=True, help="YAML aggregation manifest.")
    parser.add_argument("--artifacts-dir", required=True,
                        help="Directory written by locate_scRNA_artifacts --output-dir.")
    parser.add_argument("--analysis-id", required=True, help="Analysis identifier.")
    parser.add_argument("--reference", help="Reference name.  Skips the check that all libraries have the same one.")
    parser.add_argument("--mmc-model", help="MapMyCells model.  Default: the first model of each library.")
    parser.add_argument("--output-dir", "-o", required=True,
                        help=f"Directory for {LIBRARIES_FILE} and the join_and_filter_tsv argument files.")
    return parser.parse_args(args)


def main(args=None):
    options = parse_args(args)
    try:
        resolved = resolve(load_manifest(options.manifest), options.artifacts_dir, options.analysis_id,
                           options.reference, options.mmc_model)
    except ValueError as e:
        logger.error(str(e))
        return 1
    write_resolved(resolved, options.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
