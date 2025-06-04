#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path
import re


def count_rrna_types(gff_path):
    """Count occurrences of 5S, 16S, and 23S rRNA genes in a .rna.gff file."""
    counts = {"5S": 0, "16S": 0, "23S": 0}
    try:
        with open(gff_path, "r") as f:
            lines = [line for line in f if not line.startswith("#")]
            if not lines:
                return counts  # File is empty or has only comments
            for line in lines:
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                attributes = fields[8]
                match = re.search(r"product=([^;]+)", attributes)
                if match:
                    product = match.group(1).strip()
                    for rna_type in counts:
                        if product.startswith(rna_type):
                            counts[rna_type] += 1
    except Exception as e:
        print(f"Warning: Failed to process {gff_path}: {e}")
    return counts


def main(summary_file, result_dir, output_file):
    result_dir = Path(result_dir)
    quast_summary_df = pd.read_csv(summary_file, sep="\t", index_col=0)

    # Find all .rna.gff files and extract sample name from filename
    rna_gffs = {
        gff_path.parts[-3]: gff_path
        for gff_path in result_dir.rglob("*.rna.gff")
    }

    # Create DataFrame to hold rRNA counts
    rna_df = pd.DataFrame(columns=["5S_rRNA", "16S_rRNA", "23S_rRNA"])

    for sample_name, gff_path in rna_gffs.items():
        counts = count_rrna_types(gff_path)
        rna_df.loc[sample_name] = [counts["5S"], counts["16S"], counts["23S"]]


    # Merge and output

    merged_df = pd.concat([quast_summary_df, rna_df], axis=1)
    merged_df.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge QUAST summaries and rRNA counts.")
    parser.add_argument("--summaries", required=True, help="The summary QUAST result file.")
    parser.add_argument("--result_dir", required=True, help="Path to QUAST result directory (recursively searched).")
    parser.add_argument("--output_file", required=True, help="Path to save merged summary.")
    args = parser.parse_args()
    main(args.summaries, args.result_dir, args.output_file)
