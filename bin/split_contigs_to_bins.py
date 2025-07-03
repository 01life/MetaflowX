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

def split_contigs_to_bins(contig_fasta, contig_to_bin_file, output_dir, prefix, min_bin_length):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse contig to bin mapping
    bin_to_contigs = parse_contig_to_bin(contig_to_bin_file)
    
    # Read contig FASTA file
    contig_records = {record.id: record for record in SeqIO.parse(contig_fasta, "fasta")}
    
    # Create a FASTA file for each bin that meets the minimum length requirement
    for bin_id, contig_ids in bin_to_contigs.items():
        # Calculate total length of contigs in the bin
        total_length = sum(len(contig_records[contig_id].seq) for contig_id in contig_ids if contig_id in contig_records)
        
        if total_length >= min_bin_length:
            output_file = os.path.join(output_dir, f"{prefix}_{bin_id}.fa")
            with open(output_file, 'w') as f:
                for contig_id in contig_ids:
                    if contig_id in contig_records:
                        SeqIO.write(contig_records[contig_id], f, "fasta")
                    else:
                        print(f"Warning: Contig {contig_id} not found in FASTA file")
        else:
            print(f"Skipping bin {bin_id}: total length {total_length} bp is less than minimum {min_bin_length} bp")

def main():
    parser = argparse.ArgumentParser(description="Split contigs into bins based on contig-to-bin mapping with minimum bin length filter")
    parser.add_argument("-c", "--contig_fasta", required=True, help="Input contig FASTA file")
    parser.add_argument("-b", "--contig_to_bin", required=True, help="Contig to bin mapping TSV file")
    parser.add_argument("-o", "--output_dir", default="bins", help="Output directory for bin FASTA files (default: bins)")
    parser.add_argument("-p", "--prefix", default="All", help="Prefix for output bin FASTA file names (default: All)")
    parser.add_argument("-m", "--min_bin_length", type=int, default=2000, help="Minimum total length of bin to output (default: 10000 bp)")
    args = parser.parse_args()
    
    split_contigs_to_bins(args.contig_fasta, args.contig_to_bin, args.output_dir, args.prefix, args.min_bin_length)

if __name__ == "__main__":
    main()