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
Tables for the aggregation QC reports.  Pure pandas; nothing here plots.

Donors and cell types are labels, so callers should read them as strings (donor IDs can look numeric).
"""

import pandas as pd

from dropseq.aggregation import logger
from dropseq.aggregation.donor_cell_type_report import donor_cell_type_report
from dropseq.util.storage import default_store

CORE_POOLED_COLUMN = "median_umis_per_core_cell_type"
CORE_TYPE_COLUMN_PREFIX = "median_umis_per_nucleus_"


def read_donor_list(url, store=None):
    """
    :param url: local path or gs:// URL of a headerless file with one donor per line.  Blank lines are ignored.
    :return: list of donors in file order
    """
    store = store or default_store(url)
    return [line.strip() for line in store.read_text(url).splitlines() if line.strip()]


def read_core_cell_types(file):
    """
    :param file: headerless file with one core cell type per line.  Blank lines are ignored.
    :return: list of cell types in file order, without duplicates
    """
    with open(file) as f:
        return list(dict.fromkeys(line.strip() for line in f if line.strip()))


def expected_donors_from_files(donor_files, store=None):
    """
    The donors expected in the experiment, from the donor files of its libraries.

    :param donor_files: local paths or gs:// URLs of donor files, one per library or fewer.  The files are compared by
    content, not by location.
    :return: the sorted expected donors, or None if there are no donor files or they do not list the same donors,
    which is logged
    """
    if not donor_files:
        logger.info("No donor file; skipping donor QC")
        return None
    donor_sets = {frozenset(read_donor_list(donor_file, store)) for donor_file in donor_files}
    if len(donor_sets) > 1:
        logger.info(f"The donor files list different donors ({len(donor_sets)} distinct sets); skipping donor QC")
        return None
    return sorted(next(iter(donor_sets)))


def resolve_conditions(cell_metadata, library_column="PREFIX", conditions_file=None, condition_column=None,
                       libraries=None):
    """
    Reduce the condition, a property of the library, to a library -> condition table.

    :param conditions_file: TSV with a header and exactly two columns, (library, condition).  Column names are not
    checked.  Every library of the cell metadata, and of libraries, must be in it.
    :param condition_column: cell metadata column, which must have one value per library.
    :param libraries: more libraries, for example from the library metrics, that need a condition from conditions_file.
    Not used with condition_column, which can only describe the libraries of the cell metadata.
    :return: DataFrame with columns library, condition with one row per library in the cell metadata, or None if
    neither option was given
    :raises ValueError: if both options are given, a library has no condition, or has more than one
    """
    if conditions_file is not None and condition_column is not None:
        raise ValueError("Only one of conditions_file and condition_column may be given")
    extra_libraries = [str(library) for library in (libraries if libraries is not None else [])]
    libraries = list(dict.fromkeys(cell_metadata[library_column].dropna().astype(str)))
    if conditions_file is not None:
        table = pd.read_csv(conditions_file, sep="\t", dtype=str)
        if table.shape[1] != 2:
            raise ValueError(f"{conditions_file} must have exactly 2 columns (library, condition), "
                             f"found {table.shape[1]}")
        table.columns = ["library", "condition"]
        table = table.dropna()
        conflicting = table.groupby("library")["condition"].nunique()
        if (conflicting > 1).any():
            raise ValueError(f"Libraries with more than one condition in {conditions_file}: "
                             f"{sorted(conflicting[conflicting > 1].index)}")
        conditions = table.drop_duplicates("library").set_index("library")["condition"]
        libraries = list(dict.fromkeys([*libraries, *extra_libraries]))
        missing = [library for library in libraries if library not in conditions.index]
        if missing:
            raise ValueError(f"Libraries missing from {conditions_file}: {missing}")
    elif condition_column is not None:
        if condition_column not in cell_metadata.columns:
            raise ValueError(f"No {condition_column} column in the cell metadata")
        cells = cell_metadata.assign(**{library_column: cell_metadata[library_column].astype(str),
                                        condition_column: cell_metadata[condition_column].astype("string")})
        values = cells.dropna(subset=[condition_column]).groupby(library_column)[condition_column].unique()
        varying = sorted(library for library, v in values.items() if len(v) > 1)
        if varying:
            raise ValueError(f"{condition_column} has more than one value within librar"
                             f"{'y' if len(varying) == 1 else 'ies'} {varying}")
        missing = [library for library in libraries if library not in values.index]
        if missing:
            raise ValueError(f"No {condition_column} value for librar{'y' if len(missing) == 1 else 'ies'} {missing}")
        conditions = values.map(lambda v: v[0])
    else:
        return None
    return pd.DataFrame({"library": libraries, "condition": [conditions[library] for library in libraries]})


def donor_qc_table(cell_metadata, expected_donors, core_cell_types=(), cell_type_column=None, donor_column="donor",
                   umi_column="num_retained_transcripts"):
    """
    One row per donor, sorted by donor: the expected donors plus any donor found in the data.

    :param expected_donors: donors of the experiment.  Those with no cells get 0 nuclei and NA values.
    :param core_cell_types: cell types whose per-donor median UMIs are reported, pooled (CORE_POOLED_COLUMN) and for
    each one present in the data (CORE_TYPE_COLUMN_PREFIX + type).  Types that are absent are logged.  Needs
    cell_type_column.
    :return: DataFrame with columns donor, num_nuclei, fraction_nuclei (of nuclei with a donor), total_umis_m,
    median_umis_per_nucleus, the core cell type columns, and expected.
    """
    cells = cell_metadata.dropna(subset=[donor_column])
    total = len(cells)
    present_core = []
    if core_cell_types:
        present = set(cells[cell_type_column].dropna()) if cell_type_column else set()
        present_core = [t for t in core_cell_types if t in present]
        absent = [t for t in core_cell_types if t not in present]
        if absent:
            logger.info(f"Core cell types not found in the data: {absent}")
    donors = sorted(set(expected_donors) | set(cells[donor_column]))
    expected = set(expected_donors)
    by_donor = {donor: donor_cells for donor, donor_cells in cells.groupby(donor_column)}
    rows = []
    for donor in donors:
        donor_cells = by_donor.get(donor)
        row = {"donor": donor, "num_nuclei": 0 if donor_cells is None else len(donor_cells)}
        if donor_cells is None:
            row.update(fraction_nuclei=0.0, total_umis_m=float("nan"), median_umis_per_nucleus=float("nan"))
        else:
            umis = donor_cells[umi_column]
            row.update(fraction_nuclei=round(len(donor_cells) / total, 3), total_umis_m=round(umis.sum() / 1e6, 3),
                       median_umis_per_nucleus=umis.median())
        if present_core:
            if donor_cells is None:
                row.update({c: float("nan") for c in
                            [CORE_POOLED_COLUMN, *[CORE_TYPE_COLUMN_PREFIX + t for t in present_core]]})
            else:
                labels = donor_cells[cell_type_column]
                row[CORE_POOLED_COLUMN] = donor_cells.loc[labels.isin(present_core), umi_column].median()
                for core_type in present_core:
                    row[CORE_TYPE_COLUMN_PREFIX + core_type] = donor_cells.loc[labels == core_type,
                                                                              umi_column].median()
        row["expected"] = donor in expected
        rows.append(row)
    return pd.DataFrame(rows)


def library_cell_summary(cell_metadata, cell_type_column=None, library_column="PREFIX",
                         umi_column="num_retained_transcripts"):
    """
    Per library, after aggregation filtering.

    :return: DataFrame with columns library, num_nuclei, median_umis_per_nucleus and, if cell_type_column is given,
    num_labelled_nuclei and median_umis_per_labelled_nucleus (NaN if no nucleus has a label), in order of first
    appearance of the library
    """
    rows = []
    for library, cells in cell_metadata.groupby(library_column, sort=False):
        row = {"library": str(library), "num_nuclei": len(cells), "median_umis_per_nucleus": cells[umi_column].median()}
        if cell_type_column is not None:
            labelled = cells.dropna(subset=[cell_type_column])
            row["num_labelled_nuclei"] = len(labelled)
            row["median_umis_per_labelled_nucleus"] = labelled[umi_column].median() if len(labelled) else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def complete_donor_cell_type(cell_metadata, cell_type_column, expected_donors, donor_column="donor",
                             umi_column="num_retained_transcripts"):
    """
    The donor_cell_type_report of the cell metadata, with a row for every donor x cell type.

    :return: DataFrame with the columns of donor_cell_type_report plus expected.  Combinations with no nuclei have
    num_nuclei 0 and fraction_nuclei 0, and NA median_umis_per_nucleus (a median of nothing is not 0).  Donors are the
    expected ones plus those found; both donors and cell types are sorted.
    """
    report = donor_cell_type_report(cell_metadata, cell_type_column, donor_column, umi_column)
    donors = sorted(set(expected_donors) | set(cell_metadata[donor_column].dropna()))
    cell_types = sorted(report["cell_type"].unique())
    grid = pd.MultiIndex.from_product([donors, cell_types], names=["donor", "cell_type"])
    result = report.set_index(["donor", "cell_type"]).reindex(grid)
    result[["num_nuclei", "fraction_nuclei"]] = result[["num_nuclei", "fraction_nuclei"]].fillna(0)
    result["num_nuclei"] = result["num_nuclei"].astype(int)
    result = result.reset_index()
    result["expected"] = result["donor"].isin(set(expected_donors))
    return result


# Column names of the cell metadata that the reports work on; see normalize_cell_metadata.
LIBRARY, DONOR, CELL_TYPE, UMIS, CONDITION = "library", "donor", "cell_type", "umis", "condition"


def normalize_cell_metadata(cell_metadata, library_column="PREFIX", donor_column="donor", cell_type_column=None,
                            umi_column="num_retained_transcripts", conditions=None):
    """
    :param conditions: library -> condition table from resolve_conditions, or None
    :return: DataFrame with one row per cell and the columns LIBRARY, UMIS and, if they are available, DONOR,
    CELL_TYPE and CONDITION.  Library, donor and cell type are strings.
    """
    columns = {LIBRARY: cell_metadata[library_column].astype(str), UMIS: cell_metadata[umi_column]}
    if donor_column in cell_metadata.columns:
        columns[DONOR] = cell_metadata[donor_column].astype("string").astype(object)
    if cell_type_column is not None:
        columns[CELL_TYPE] = cell_metadata[cell_type_column].astype("string").astype(object)
    cells = pd.DataFrame(columns)
    if conditions is not None:
        cells[CONDITION] = cells[LIBRARY].map(dict(zip(conditions["library"], conditions["condition"])))
    return cells.reset_index(drop=True)


def cell_type_table(cells, group_columns=()):
    """
    The number and fraction of the nuclei of each group that are of each cell type, and their median UMIs.

    :param cells: normalized cell metadata with CELL_TYPE.  Cells with a missing group value are ignored.
    :param group_columns: columns that define the groups, for example [LIBRARY] or [CONDITION].  Empty for all cells.
    :return: DataFrame with the group columns, CELL_TYPE, num_nuclei, fraction_nuclei (of all nuclei of the group,
    including those with no cell type) and median_umis_per_nucleus, for every group x cell type.  A group with no
    nuclei of a cell type has num_nuclei 0, fraction_nuclei 0 and an NA median.
    """
    keys = list(group_columns) or ["_all"]
    cells = cells.dropna(subset=list(group_columns))
    if not group_columns:
        cells = cells.assign(_all="all")
    totals = cells.groupby(keys).size()
    typed = cells.dropna(subset=[CELL_TYPE])
    stats = typed.groupby(keys + [CELL_TYPE]).agg(num_nuclei=(CELL_TYPE, "size"),
                                                  median_umis_per_nucleus=(UMIS, "median"))
    cell_types = sorted(typed[CELL_TYPE].unique())
    groups = [g if isinstance(g, tuple) else (g,) for g in totals.index]
    grid = pd.MultiIndex.from_tuples([(*g, t) for g in groups for t in cell_types], names=keys + [CELL_TYPE])
    result = stats.reindex(grid)
    result["num_nuclei"] = result["num_nuclei"].fillna(0).astype(int)
    result["fraction_nuclei"] = result["num_nuclei"] / totals.reindex(grid.droplevel(CELL_TYPE)).to_numpy()
    result = result.reset_index()
    return result.drop(columns="_all") if not group_columns else result


def donor_library_counts(cells):
    """
    :return: DataFrame with LIBRARY, DONOR, num_nuclei (and CONDITION if there is one) of the cells with a donor
    """
    keys = [LIBRARY, DONOR] if CONDITION not in cells.columns else [LIBRARY, CONDITION, DONOR]
    return cells.dropna(subset=[DONOR]).groupby(keys).size().rename("num_nuclei").reset_index()
