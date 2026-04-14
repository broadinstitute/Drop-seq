import argparse
import re
import sys
from typing import TextIO, Optional

import pandas as pd


# Write Sample IDs and optionally Sample Names from a barcodes file
# The Sample IDs are derived from the SampleName column by replacing
# non-alphanumeric characters with underscores.

def get_args():
    parser = argparse.ArgumentParser(description="Write the list of sample IDs from the barcodes file.")
    parser.add_argument(
        "-b",
        "--barcodes",
        required=True,
        help="Tab separated file with a header and at least a column named SampleName",
    )
    parser.add_argument(
        "--ids",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="Output file for the list of sample ids.  Default: stdout",
    )
    parser.add_argument(
        "--names",
        help="Optional output file for the list of sample names.",
    )

    args = parser.parse_args()

    return args


def write_sample_ids(
        barcodes: pd.DataFrame,
        ids_output: TextIO,
        names_file: Optional[str],
) -> None:
    """
    Write the sample ids and names.
    """
    assert "SampleName" in barcodes.columns

    if names_file:
        with open(names_file, "w") as names_output:
            for sample_name in barcodes.SampleName:
                names_output.write(f"{sample_name}\n")
    for sample_name in barcodes.SampleName:
        sample_id = re.sub(r"[^a-zA-Z0-9_-]", "_", sample_name)
        ids_output.write(f"{sample_id}\n")


def main() -> None:
    args = get_args()

    barcodes_df = pd.read_csv(args.barcodes, sep="\t")
    write_sample_ids(barcodes_df, args.ids, args.names)


if __name__ == "__main__":
    main()
