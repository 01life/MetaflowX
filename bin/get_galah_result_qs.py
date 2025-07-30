#!/usr/bin/env python

import argparse
import pandas as pd
import os

def extract_bin_ids_from_fa(input_dir):
    """
    Recursively find all .fa files in the input directory and return a set of filenames without the extension.
    """
    bin_ids = set()
    for root, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith(".fa"):
                bin_ids.add(os.path.splitext(f)[0])
    return bin_ids

def filter_checkm2_results(checkm2_file, bin_ids):
    """
    Read the checkm2 TSV result file and filter rows where the 'Name' field ends with any of the bin IDs.
    """
    df = pd.read_csv(checkm2_file, sep='\t')
    df_filtered = df[df['Name'].apply(lambda name: any(name.endswith(bin_id) for bin_id in bin_ids))]
    return df_filtered

def main():
    parser = argparse.ArgumentParser(description="Extract checkm2 results for bins found in .fa files under a directory")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing .fa files")
    parser.add_argument("-c", "--checkm2", required=True, help="Path to checkm2 result TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to output filtered result TSV file")

    args = parser.parse_args()

    # Extract bin IDs from .fa filenames
    bin_ids = extract_bin_ids_from_fa(args.input_dir)
    print(f"[INFO] Found {len(bin_ids)} bin(s) in {args.input_dir}: {sorted(bin_ids)}")

    # Filter checkm2 results based on extracted bin IDs
    filtered_df = filter_checkm2_results(args.checkm2, bin_ids)
    print(f"[INFO] Matched {len(filtered_df)} entries, saving to {args.output}")

    # Save filtered results to output file
    filtered_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
