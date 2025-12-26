import argparse
import os
import re
import sys
from typing import Optional

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Subset fastq files for a specific sample.")
    indexes_group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument(
        "-f",
        "--fastqs",
        required=True,
        help="File containing the paths to the fastq files.  This file should contain one path per line.",
    )
    parser.add_argument(
        "-s",
        "--sample-name",
        required=True,
        help="sample_name to filter the barcodes to.",
    )
    indexes_group.add_argument(
        "-n",
        "--sample-number",
        type=int,
        help="Sample number to filter the barcodes to.  This is the sample number in the barcodes file, starting at 1.",
    )
    indexes_group.add_argument(
        "-b",
        "--barcodes",
        help="Tab separated file with a header and at least a column named SampleName",
    )
    parser.add_argument(
        "-l",
        "--lanes",
        dest="lanes_csv",
        default="1,2,3,4,5,6,7,8",
        help="Comma-separated list of lanes to include in the file of file names.  Default: %(default)s",
    )
    parser.add_argument(
        "-r",
        "--reads",
        dest="reads_csv",
        default="1,2",
        help="Comma-separated list of reads to include in the file of file names.  Default: %(default)s",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="Output file.  Default: stdout",
    )

    args = parser.parse_args()

    lanes_split = args.lanes_csv.split(",")

    if len(lanes_split) > 8:
        parser.error("--lanes must contain 8 or fewer comma-separated int values")

    try:
        args.lanes = [int(lane) for lane in lanes_split]
    except ValueError:
        parser.error("--lanes must contain 8 or fewer comma-separated int values")

    for lane in args.lanes:
        if lane < 1 or lane > 8:
            parser.error("--lanes must contain int values between 1 and 8")

    reads_split = args.reads_csv.split(",")

    if len(reads_split) > 2:
        parser.error("--reads must contain 2 or fewer comma-separated int values")

    try:
        args.reads = [int(read) for read in reads_split]
    except ValueError:
        parser.error("--reads must contain 2 or fewer comma-separated int values")

    for read in args.reads:
        if read < 1 or read > 2:
            parser.error("--reads must contain int values of 1 and/or 2")

    return args


def list_sample_fastqs(
        fastq_paths: dict[str, str],
        lanes: list[int],
        reads: list[int],
        sample_id: str,
        sample_number: int,
) -> list[str]:
    sample_fastqs = []
    for lane in lanes:
        for read in reads:
            file_name = f"{sample_id}_S{sample_number}_L00{lane}_R{read}_001.fastq.gz"
            if file_name not in fastq_paths:
                raise FileNotFoundError(f"File {file_name} not found.")
            sample_fastq = fastq_paths[file_name]
            sample_fastqs.append(sample_fastq)
    return sample_fastqs


def subset_sample_fastqs(
        fastqs: list[str],
        selected_sample_name: str,
        selected_sample_number: Optional[int],
        lanes: list[int],
        reads: list[int],
        barcodes: Optional[pd.DataFrame],
) -> list[str]:
    """
    Write the file names for the specified lanes and reads.
    """
    fastq_paths = {}
    for fastq in fastqs:
        filename = os.path.basename(fastq)
        fastq_paths[filename] = fastq

    if barcodes is not None:
        assert "SampleName" in barcodes.columns
        for sample_index, sample_name in enumerate(barcodes.SampleName):

            sample_id = re.sub(r"[^a-zA-Z0-9_-]", "_", sample_name)

            if sample_name != selected_sample_name and sample_id != selected_sample_name:
                continue

            return list_sample_fastqs(fastq_paths, lanes, reads, sample_id, sample_index + 1)
        return []
    elif selected_sample_number:
        selected_sample_id = re.sub(r"[^a-zA-Z0-9_-]", "_", selected_sample_name)
        return list_sample_fastqs(fastq_paths, lanes, reads, selected_sample_id, selected_sample_number)
    else:
        raise ValueError("Either barcodes or sample_number must be provided.")


def main() -> None:
    args = get_args()

    barcodes_df = None
    if args.barcodes:
        barcodes_df = pd.read_csv(args.barcodes, sep="\t")

    with open(args.fastqs, "r") as f:
        fastqs = [line.strip() for line in f if line.strip()]

    sample_fastqs = subset_sample_fastqs(
        fastqs,
        args.sample_name,
        args.sample_number,
        args.lanes,
        args.reads,
        barcodes_df,
    )

    if not sample_fastqs:
        print(f"No fastq files found for sample {args.sample_name}", file=sys.stderr)
        sys.exit(1)

    for sample_fastq in sample_fastqs:
        args.output.write(f"{sample_fastq}\n")


if __name__ == "__main__":
    main()
