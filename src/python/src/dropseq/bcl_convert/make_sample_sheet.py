import argparse
import os

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Generate sample sheet to convert bcl files to fastq files")
    parser.add_argument(
        "-b",
        "--barcodes",
        required=True,
        help="Tab separated file with a header of columns SampleName, IndexBarcode1, IndexBarcode2",
    )
    parser.add_argument(
        "-s",
        "--sample-sheet",
        required=True,
        help="Path to output the sample sheet file",
    )
    parser.add_argument(
        "--rcib1",
        default=False,
        action="store_true",
        help="Reverse complement the IndexBarcode1 from the barcodes file to create the index1 on the sample sheet",
    )
    parser.add_argument(
        "--rcib2",
        default=False,
        action="store_true",
        help="Reverse complement the IndexBarcode2 from the barcodes file to create the index2 on the sample sheet",
    )

    args = parser.parse_args()

    return args


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a sequence.
    """
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }
    return "".join(complement[base] for base in reversed(seq))


def barcodes_to_sample_sheet(
        barcodes: pd.DataFrame,
        reverse_complement_index1: bool,
        reverse_complement_index2: bool,
) -> pd.DataFrame:
    """
    Read the barcode indexes from a tsv returning a dataframe with columns Sample_ID, index, and maybe index2.
    """
    assert "SampleName" in barcodes.columns
    assert "IndexBarcode1" in barcodes.columns
    if not "IndexBarcode2" in barcodes.columns:
        barcodes["IndexBarcode2"] = None
        is_single_index = True
    else:
        is_single_index = barcodes["IndexBarcode2"].isna().all()
    barcodes.rename(
        columns={"SampleName": "Sample_ID", "IndexBarcode1": "index", "IndexBarcode2": "index2"},
        inplace=True,
    )
    # In 'SampleName' replace anything not alphanumeric, -, or _ with '_'
    barcodes["Sample_ID"] = barcodes["Sample_ID"].str.replace(r"[^a-zA-Z0-9-_]", "_", regex=True)

    if reverse_complement_index1:
        barcodes["index"] = barcodes["index"].apply(reverse_complement)

    if not is_single_index and reverse_complement_index2:
        barcodes["index2"] = barcodes["index2"].apply(reverse_complement)

    # Reorder columns
    barcodes["Sample_Project"] = barcodes["Sample_ID"]
    if is_single_index:
        barcodes = barcodes[["Sample_Project", "Sample_ID", "index"]]
    else:
        barcodes = barcodes[["Sample_Project", "Sample_ID", "index", "index2"]]

    return barcodes


def write_sample_sheet(
        sample_sheet_df: pd.DataFrame,
        sample_sheet_path: str,
) -> None:
    """
    Write the sample sheet to a file.
    """
    sheet_path_dir = os.path.dirname(sample_sheet_path)
    if sheet_path_dir:
        os.makedirs(sheet_path_dir, exist_ok=True)
    with open(sample_sheet_path, 'w') as f:
        f.write("[Header]\n")
        f.write("FileFormatVersion,2\n")
        f.write("\n")
        f.write("[BCLConvert_Settings]\n")
        f.write("\n")
        f.write("[BCLConvert_Data]\n")
    sample_sheet_df.to_csv(sample_sheet_path, mode='a', index=False)


def main() -> None:
    args = get_args()

    print(f"Barcodes: {args.barcodes}")
    print(f"Reverse Complement IndexBarcode1: {args.rcib1}")
    print(f"Reverse Complement IndexBarcode2: {args.rcib2}")
    print(f"Sample sheet: {args.sample_sheet}")

    # Get the barcode indexes from barcodes tsv
    barcodes_df = pd.read_csv(args.barcodes, sep="\t")
    sample_sheet_df = barcodes_to_sample_sheet(barcodes_df, args.rcib1, args.rcib2)
    write_sample_sheet(sample_sheet_df, args.sample_sheet)


if __name__ == "__main__":
    main()
