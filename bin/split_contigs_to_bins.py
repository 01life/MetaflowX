#!/usr/bin/env python

import argparse
from Bio import SeqIO
import os

def parse_contig_to_bin(file_path):
    bin_to_contigs = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            bin_id, contig_id = line.strip().split('\t')
            if bin_id not in bin_to_contigs:
                bin_to_contigs[bin_id] = []
            bin_to_contigs[bin_id].append(contig_id)
    return bin_to_contigs

def split_contigs_to_bins(contig_fasta, contig_to_bin_file, output_dir, prefix):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse contig to bin mapping
    bin_to_contigs = parse_contig_to_bin(contig_to_bin_file)
    
    # Read contig FASTA file
    contig_records = {record.id: record for record in SeqIO.parse(contig_fasta, "fasta")}
    
    # Create a FASTA file for each bin
    for bin_id, contig_ids in bin_to_contigs.items():
        output_file = os.path.join(output_dir, f"{prefix}_{bin_id}.fa")
        with open(output_file, 'w') as f:
            for contig_id in contig_ids:
                if contig_id in contig_records:
                    SeqIO.write(contig_records[contig_id], f, "fasta")
                else:
                    print(f"Warning: Contig {contig_id} not found in FASTA file")

def main():
    parser = argparse.ArgumentParser(description="Split contigs into bins based on contig-to-bin mapping")
    parser.add_argument("-c", "--contig_fasta", required=True, help="Input contig FASTA file")
    parser.add_argument("-b", "--contig_to_bin", required=True, help="Contig to bin mapping TSV file")
    parser.add_argument("-o", "--output_dir", default="bins", help="Output directory for bin FASTA files (default: bins)")
    parser.add_argument("-p", "--prefix", default="All", help="Prefix for output bin FASTA file names (default: All)")
    args = parser.parse_args()
    
    split_contigs_to_bins(args.contig_fasta, args.contig_to_bin, args.output_dir, args.prefix)

if __name__ == "__main__":
    main()