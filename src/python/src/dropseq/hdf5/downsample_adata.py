#!/usr/bin/env python3
# MIT License
#
# Copyright 2024 Broad Institute
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
from tqdm import tqdm
import anndata as ad
import pandas as pd
from typing import List, Optional
import argparse
import os
import sys
import numpy as np

import dropseq.hdf5.filters as filters

import logging
logger = logging.getLogger(__name__)

def compute_downsampling_counts(
    h5ad_files: List[str],
    group_cols: List[str],
    max_cells: int,
    expressions: List[str]
) -> pd.DataFrame:
    """
    Compute the exact number of cells to retain for each group in each file,
    such that the total retained per group across all files does not exceed `max_cells`.

    Fully filtered files will still be represented in the output with "NA" values in group_cols
    and all counts set to 0.

    Parameters
    ----------
    h5ad_files : str
        List of .h5ad files.
    group_cols : List[str]
        Column names in .obs to use for grouping.
    max_cells : int
        Maximum number of cells to retain per group across all files.
    expressions : List[str]
        Python expressions to filter cells before grouping.

    Returns
    -------
    pd.DataFrame
        A table with one row per file/group combination, including:
            - group_cols
            - 'filename': the file this row is from
            - 'count': number of filtered cells observed for this group in the file
            - 'total': total observed for the group across all files
            - 'target_total': capped target for group total (min(total, max_cells))
            - 'expected_final_count': exact number of cells to retain from this file
    """
    all_counts = []
    total_cells = 0

    for path in tqdm(h5ad_files, desc="Building Downsampling Plan"):
        fname=os.path.basename(path)
        adata = ad.read_h5ad(path, backed="r")
        adata = filters.evaluate_filters(adata, expressions)
        n_total = adata.n_obs
        total_cells += n_total

        # If no cells pass filters, emit a zero-row placeholder
        if adata.n_obs == 0:
            empty_row = {
                **{col: "NA" for col in group_cols},
                "filename": fname,
                "count": 0
            }
            all_counts.append(pd.DataFrame([empty_row]))
            continue

        df = adata.obs[group_cols].copy()
        df[group_cols] = df[group_cols].astype(str).replace("nan", "NA")

        group_sizes = df.groupby(group_cols, observed=True).size().reset_index(name='count')
        group_sizes['filename'] = fname
        all_counts.append(group_sizes)

    per_file_counts = pd.concat(all_counts)

    # Total observed cells per group across all files
    global_totals = (
        per_file_counts
        .groupby(group_cols, observed=True)["count"]
        .sum()
        .reset_index(name="total")
    )

    per_file_counts = per_file_counts.merge(global_totals, on=group_cols, how="left")
    per_file_counts["total"] = per_file_counts["total"].fillna(0).astype(int)

    # Cap the total per group at max_cells (e.g., 1000)
    per_file_counts["target_total"] = per_file_counts["total"].clip(upper=max_cells)

    # Avoid division by zero
    safe_total = per_file_counts["total"].replace(0, np.nan)

    # Compute proportional downsampling target per file, round up to nearest integer.
    # adjust_expected_counts will fix any integer overflow later.
    per_file_counts["expected_final_count"] = (
            per_file_counts["count"] * (per_file_counts["target_total"] / safe_total)
    ).fillna(0).apply(np.ceil).clip(upper=per_file_counts["count"]).astype(int)

    per_file_counts= adjust_expected_counts(per_file_counts, group_cols, max_cells)

    logger.info(
        "Saw %d cells across all files. Expected final total after downsampling: %d",
        total_cells,
        per_file_counts["expected_final_count"].sum()
    )

    return per_file_counts

def adjust_expected_counts(df: pd.DataFrame, group_cols: List[str], max_cells: int) -> pd.DataFrame:
    """
    Adjusts expected_final_count so that each group's sum does not exceed max_cells.
    """
    adjusted = df.copy()
    grouped = adjusted.groupby(group_cols, observed=True)

    for name, group in grouped:
        total_expected = group['expected_final_count'].sum()

        if total_expected > max_cells:
            excess = total_expected - max_cells

            # Sort rows with room to decrement (expected_final_count > 1)
            candidates = group[group["expected_final_count"] > 1].copy()
            candidates = candidates.sort_values("expected_final_count", ascending=False)

            # Use index to track updates in adjusted
            decrement_pool = candidates.index.tolist()
            i = 0

            while excess > 0:
                idx = decrement_pool[i % len(decrement_pool)]
                if adjusted.at[idx, "expected_final_count"] > 1:
                    adjusted.at[idx, "expected_final_count"] -= 1
                    excess -= 1
                i += 1

    return adjusted




def downsample_anndata(
        adata: ad.AnnData,
        group_cols: List[str],
        count_df: pd.DataFrame,
        filename: str,
        seed: Optional[int] = None
) -> ad.AnnData:
    """
    Downsample an AnnData object to retain exactly expected_final_count cells per group,
    based on a precomputed count plan.
    """
    rng = np.random.default_rng(seed)

    obs = adata.obs[group_cols].copy()
    obs[group_cols] = obs[group_cols].astype(str).replace("nan", "NA")

    # Preserve original index for subsetting later
    obs = obs.reset_index(drop=False)

    plan = count_df[count_df["filename"] == filename]
    if plan.empty:
        raise ValueError(f"No downsampling plan found for file: {filename}")

    plan = plan[group_cols + ['expected_final_count']]

    merged = obs.merge(plan, on=group_cols, how="left", indicator=True)

    unmatched = merged['_merge'] == 'left_only'
    if unmatched.any():
        bad = merged.loc[unmatched, group_cols].drop_duplicates()
        raise ValueError(
            f"The following groups in the AnnData object were not found in the count plan for {filename}:\n"
            f"{bad.to_string(index=False)}"
        )

    keep_indices = []
    for _, group in merged.groupby(group_cols, observed=True):
        n_expected = group["expected_final_count"].iloc[0]
        if n_expected >= len(group):
            keep_indices.extend(group["index"].tolist())
        else:
            keep_indices.extend(
                rng.choice(group["index"].to_numpy(), size=n_expected, replace=False)
            )

    return adata[keep_indices]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute downsampling rates for grouped cells across multiple AnnData files in a directory."
    )

    parser.add_argument(
        "--dir",
        required=True,
        help="Directory containing .h5ad files."
    )

    parser.add_argument(
        "--group-cols",
        action="append",
        required=True,
        help="Repeatable: column(s) in .obs to use for grouping (e.g. --group-cols donor_id)"
    )

    parser.add_argument(
        "--max-cells",
        type=int,
        required=True,
        help="Maximum number of cells to retain per group across all files."
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional path to write the output CSV file with downsampling rates."
    )

    parser.add_argument("--expr", action='append', default=[],
                        help="Optional Python expressions to filter cells (e.g., 'adata.obs[\"n_genes\"] > 500').")

    # Show help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def main():
    args = parse_args()

    if not any("ipykernel" in m for m in sys.modules):
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s] %(message)s",
                level=logging.INFO,
                datefmt='%Y-%m-%d %H:%M:%S'
            )

    rate_df = compute_downsampling_counts(
        input_dir=args.dir,
        group_cols=args.group_cols,
        max_cells=args.max_cells,
        expressions=args.expr
    )

    if args.output:
        rate_df.to_csv(args.output, index=False, header=True, sep="\t")
        print(f"Saved downsampling rates to {args.output}")
    else:
        print(rate_df)

if __name__ == "__main__":
    main()
