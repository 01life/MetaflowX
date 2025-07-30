#!/usr/bin/env python


import pandas as pd
import argparse
import os
from pathlib import Path
from typing import Generator

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge multiple Kraken2 result files into a single table.")
    parser.add_argument(
        "-i", "--input-files",
        nargs='+',
        required=True,
        help="List of Kraken2 result files to merge (tab-separated). Supports wildcards (e.g., data/*.tsv)."
    )
    parser.add_argument(
        "-o", "--output-file",
        required=True,
        help="Output file path for the merged table."
    )
    return parser.parse_args()

def read_kraken_file(filepath: str) -> pd.DataFrame:
    """Read a Kraken2 result file and return a DataFrame with clade_name as index and sample column."""
    # Use only necessary columns, assume count column is the second
    df = pd.read_csv(filepath, sep='\t', usecols=[0, 1])
    
    # # Rename the count column to sample ID (assume filename without extension)
    # sample_id = Path(filepath).stem
    # df.columns = ['clade_name', sample_id]
    
    return df.set_index('clade_name')

def merge_kraken_files(filepaths: list[str]) -> pd.DataFrame:
    """Merge Kraken2 files one-by-one to control memory usage."""
    merged_df = None

    for i, path in enumerate(filepaths):
        print(f"[{i+1}/{len(filepaths)}] Processing: {path}")
        df = read_kraken_file(path)

        if merged_df is None:
            merged_df = df
        else:
            merged_df = merged_df.join(df, how='outer')
    
    # Replace NaN with 0
    merged_df = merged_df.fillna(0)

    # Convert to int if needed
    for col in merged_df.columns:
        merged_df[col] = merged_df[col].astype(int)

    # Reset index
    return merged_df.reset_index().sort_values("clade_name")

def main():
    args = parse_arguments()
    
    # Expand wildcards
    expanded_files = []
    for pattern in args.input_files:
        expanded = list(Path('.').glob(pattern))
        if not expanded:
            raise FileNotFoundError(f"No files matched pattern: {pattern}")
        expanded_files.extend(expanded)

    if not expanded_files:
        raise ValueError("No valid input files provided.")

    # Ensure output directory exists
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Merge and write output
    merged_df = merge_kraken_files([str(p) for p in expanded_files])
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"✅ Merged result saved to: {output_path}")

if __name__ == "__main__":
    main()
