#!/usr/bin/env python

import pandas as pd
import argparse

def read_pfam_tblout(filepath: str) -> pd.DataFrame:
    """Read hmmsearch --tblout from Pfam and extract $1, $4, $5"""
    records = []
    try:
        with open(filepath) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 6:
                    records.append((parts[0], parts[3], parts[4]))  # $1 = target, $4 = query, $5 = score
    except FileNotFoundError:
        print(f"[Warning] Pfam file not found: {filepath}")
    return pd.DataFrame(records, columns=["target", "query", "score"])

def read_tigr_tblout(filepath: str) -> pd.DataFrame:
    """Read hmmsearch --tblout from TIGR and extract $1, $3, $5"""
    records = []
    try:
        with open(filepath) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 6:
                    records.append((parts[0], parts[2], parts[4]))  # $1 = target, $3 = query, $5 = score
    except FileNotFoundError:
        print(f"[Warning] TIGR file not found: {filepath}")
    return pd.DataFrame(records, columns=["target", "query", "score"])

def main():
    parser = argparse.ArgumentParser(description="Merge Pfam and TIGR hmmsearch tblout results")
    parser.add_argument("--pfam", required=True, help="Path to Pfam hmmsearch tblout file")
    parser.add_argument("--tigr", required=True, help="Path to TIGR hmmsearch tblout file")
    parser.add_argument("--out", required=True, help="Path to output merged file")
    args = parser.parse_args()

    df_pfam = read_pfam_tblout(args.pfam)
    df_tigr = read_tigr_tblout(args.tigr)

    merged_df = pd.concat([df_pfam, df_tigr], ignore_index=True)

    merged_df.to_csv(args.out, sep="\t", header=False, index=False)
    print(f"[Done] Merged results written to: {args.out} ({len(merged_df)} hits)")

if __name__ == "__main__":
    main()
