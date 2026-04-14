import argparse

import pandas as pd


# Find the best orientation for index barcodes based on multiple % Undetermined reads results.
#
# Works around any issues with bcl-convert, RunInfo.xml, isReverseComplement, and barcode orientations which
# "can be wrong in some cases due to software bugs."[0]
#
# See:
#  - [0] https://help.dragen.illumina.com/product-guide/dragen-v4.4/bcl-conversion#barcode-index-orientation-for-index-and-index2-columns
#  - [1] https://knowledge.illumina.com/software/general/software-general-faq-list/000005078


def get_args():
    parser = argparse.ArgumentParser(description="Find the best barcode orientation based on % Undetermined reads.")
    parser.add_argument(
        "-d",
        "--demux",
        required=True,
        help="Comma-separated list of Demultiplex_Stats.csv files"
    )
    parser.add_argument(
        "--rcib1",
        required=True,
        help="Comma-separated list of booleans indicating whether to reverse complement index barcode 1 for each lane",
    )
    parser.add_argument(
        "--rcib2",
        help="Comma-separated list of booleans indicating whether to reverse complement index barcode 2 for each lane",
    )
    parser.add_argument(
        "--best-rcib1",
        required=True,
        help="Output path for the index barcode 1 orientation result",
    )
    parser.add_argument(
        "--best-rcib2",
        help="Output path for the index barcode 2 orientation result",
    )

    args = parser.parse_args()

    args.demux_list = args.demux.split(",")
    args.rcib1_list = args.rcib1.split(",")

    if len(args.demux_list) != len(args.rcib1_list):
        parser.error("The number of --demux files must match the number of --rcib1 values")

    if args.rcib2:
        args.rcib2_list = args.rcib2.split(",")
        if len(args.demux_list) != len(args.rcib2_list):
            parser.error("The number of --demux files must match the number of --rcib2 values")
    else:
        args.rcib2_list = []

    if bool(args.rcib2) != bool(args.best_rcib2):
        parser.error("--rcib2 and --best-rcib2 must be provided together")

    return args


def find_best_barcode_orientation(
        demux_list,
        rcib1_list,
        rcib2_list,
):
    index_best = -1
    undetermined_lowest = 2.0

    for i in range(len(demux_list)):
        demux_df = pd.read_csv(demux_list[i])

        undetermined_row = demux_df[demux_df['SampleID'] == 'Undetermined']
        if not undetermined_row.empty:
            undetermined_reads = undetermined_row['% Reads'].values[0]
        else:
            undetermined_reads = 1.0 - demux_df['% Reads'].sum()

        print(
            f"Demux file: {demux_list[i]}, " +
            f"RCIB1: {rcib1_list[i]}, " +
            (f"RCIB2: {rcib2_list[i]}, " if rcib2_list else "") +
            f"% Reads Undetermined: {undetermined_reads}"
        )

        if undetermined_lowest > undetermined_reads:
            undetermined_lowest = undetermined_reads
            index_best = i

    return index_best


def write_best_orientations(
        index_best,
        rcib1_list,
        rcib2_list,
        best_rcib1_path,
        best_rcib2_path,
):
    with open(best_rcib1_path, "w") as f:
        f.write(f"{rcib1_list[index_best]}\n")

    if best_rcib2_path:
        with open(best_rcib2_path, "w") as f:
            f.write(f"{rcib2_list[index_best]}\n")


def main() -> None:
    args = get_args()

    index_best = find_best_barcode_orientation(
        args.demux_list,
        args.rcib1_list,
        args.rcib2_list,
    )

    write_best_orientations(
        index_best,
        args.rcib1_list,
        args.rcib2_list,
        args.best_rcib1,
        args.best_rcib2,
    )


if __name__ == "__main__":
    main()
