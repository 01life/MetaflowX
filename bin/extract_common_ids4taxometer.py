#!/usr/bin/env python3

import argparse
import pathlib
import gzip
from typing import Tuple, Dict, List
from Bio import SeqIO

def read_file_to_dict(file_path: pathlib.Path) -> Tuple[str, Dict[str, str], List[str]]:
    """Read a tab-separated file into a dict using the first column as key. Preserve header and ID order."""
    with file_path.open() as f:
        header = f.readline().strip()
        data = {}
        id_order = []
        for line in f:
            if line.strip():
                key = line.split('\t', 1)[0]
                data[key] = line.strip()
                id_order.append(key)
    return header, data, id_order

def write_filtered_file(output_path: pathlib.Path, header: str, data: Dict[str, str], ordered_keys: List[str]) -> None:
    """Write lines matching ordered_keys to output, preserving the header and order."""
    with output_path.open('w') as f:
        print(header, file=f)
        for key in ordered_keys:
            print(data[key], file=f)

def open_fasta_file(file_path: pathlib.Path):
    """Open FASTA file (gzipped or not) and return file handle."""
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")  # text mode
    return file_path.open()

def main():
    parser = argparse.ArgumentParser(description="Extract rows with common IDs (in same order) from 2 TSVs and a FASTA (gz or not).")
    parser.add_argument("file1", type=pathlib.Path, help="First input TSV file")
    parser.add_argument("file2", type=pathlib.Path, help="Second input TSV file")
    parser.add_argument("file3", type=pathlib.Path, help="FASTA file (contigs, supports .gz)")
    parser.add_argument("--out1", type=pathlib.Path, default="filtered_file1.tsv", help="Output for file1")
    parser.add_argument("--out2", type=pathlib.Path, default="filtered_file2.tsv", help="Output for file2")
    parser.add_argument("--out3", type=pathlib.Path, default="filtered_contig.fa", help="Output for FASTA")
    args = parser.parse_args()

    # Read file1 and file2 into dictionaries
    header1, data1, order1 = read_file_to_dict(args.file1)
    header2, data2, _ = read_file_to_dict(args.file2)

    # Keep only IDs common to both, and keep order from file1
    common_ids_ordered = [id_ for id_ in order1 if id_ in data2]

    # Write filtered TSV files
    write_filtered_file(args.out1, header1, data1, common_ids_ordered)
    write_filtered_file(args.out2, header2, data2, common_ids_ordered)

    # Read FASTA sequences (from .fa or .fa.gz)
    with open_fasta_file(args.file3) as handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    # Write filtered FASTA
    with args.out3.open('w') as f:
        for id_ in common_ids_ordered:
            if id_ in seq_dict:
                SeqIO.write(seq_dict[id_], f, "fasta")

    print(f"Found {len(common_ids_ordered)} common IDs (in order)")
    print(f"Output written to: {args.out1}, {args.out2}, {args.out3}")

if __name__ == "__main__":
    main()
