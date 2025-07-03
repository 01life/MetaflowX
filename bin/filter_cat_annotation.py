#!/usr/bin/env python

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter CAT contig2classification.txt by removing 'no taxid assigned' entries and retaining key columns."
    )
    parser.add_argument("-i", "--input", required=True, help="Input CAT contig2classification.txt file")
    parser.add_argument("-o", "--output1", required=True, help="Output file containing filtered results with scores")
    parser.add_argument("-O", "--output2", required=True, help="Output file containing filtered results without scores")
    return parser.parse_args()

def main():
    args = parse_args()

    # Read the CAT classification file
    df = pd.read_csv(
        args.input,
        sep="\t",
        comment="#",
        header=None,
        names=["contig", "classification", "reason", "lineage", "lineage_scores"],
        dtype=str  # Read all fields as strings
    )

    # Remove rows with unclassified contigs
    df_filtered = df[df["classification"] != "no taxid assigned"]

    # Retain only contig, lineage, and lineage_scores columns
    df_out = df_filtered[["contig", "lineage", "lineage_scores"]]

    df_out.columns = ["contigs", "predictions", "scores"]

    df_out_2 = df_filtered[["contig", "lineage"]]
    df_out_2.columns = ["contigs", "predictions"]


    # Write the filtered results to output
    df_out.to_csv(args.output1, sep="\t", index=False)
    df_out_2.to_csv(args.output2, sep="\t", index=False)
    print(f"Filtered {len(df) - len(df_out)} rows. Output written to: {args.output1}")

if __name__ == "__main__":
    main()
