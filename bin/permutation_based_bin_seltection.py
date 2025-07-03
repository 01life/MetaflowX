#!/usr/bin/env python

import os
from Bio import SeqIO
import argparse

def create_output_dir(output_dir):
    os.makedirs(output_dir, exist_ok=True)

def get_highest_qs_binid(contig_bin_dir, quality_dir):
    result = {}
    for contigid, binids in contig_bin_dir.items():
        if binids:
            highest_qs_binid = max(binids, key=lambda binid: quality_dir.get(binid, float('-inf')))
            result[contigid] = highest_qs_binid
    return result

def parse_fasta(fasta_file, contig_dir):
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_dir[record.id] = record.seq

# def main(contig_bin_file, checkm2_file, contig_fasta_file, output_file, output_dir,out_contigs2binF,min_completeness, max_contamination):
#     contig_bin_dir = {}
#     quality_dir = {}
#     quality_info_dir = {}
#     bin_contig_dir = {}
#     contig_seq_dir = {}

#     # Parse FASTA file and store sequences
#     parse_fasta(contig_fasta_file, contig_seq_dir)

#     # Create output directory if it doesn't exist
#     create_output_dir(output_dir)

#     # Process the input files
#     with open(contig_bin_file, 'r') as contig_bin_f, open(checkm2_file, 'r') as checkm2_f, open(output_file, 'w') as out_f, open(out_contigs2binF,'w') as out_contigs2binFile:

#         # Write header for the output file
#         header = checkm2_f.readline()
#         out_f.write(header)

#         # Read quality scores from checkm2 file
#         for line in checkm2_f:
#             fields = line.strip().split("\t")
#             binid = fields[0]
#             completeness = float(fields[1])
#             contamination = float(fields[2])

#             if completeness > min_completeness and contamination < max_contamination:

#                 qs = completeness - 5 * contamination
#                 quality_dir[binid] = qs
#                 quality_info_dir[binid] = fields

#         # Read contig-bin mappings from the contigBin file
#         for line in contig_bin_f:
#             contigid, binid = line.strip().split("\t")
#             contig_bin_dir.setdefault(contigid, []).append(binid)
#             bin_contig_dir.setdefault(binid, []).append(contigid)

#         # Get the highest quality bin IDs for each contig
#         print(len(quality_dir.keys()))
#         highest_qs_binids = get_highest_qs_binid(contig_bin_dir, quality_dir)
#         highest_qs_binids_QS_dir = {}

#         for one_contig, binid in highest_qs_binids.items():
#             if binid in quality_dir:
#                 highest_qs_binids_QS_dir[binid] = quality_dir[binid]

#         # Sort bins by quality score in descending order
#         sorted_bins = sorted(highest_qs_binids_QS_dir, key=highest_qs_binids_QS_dir.get, reverse=True)
#         print(sorted_bins)
#         best_bins = []
#         existing_contig_list = []

#         # Select bins based on coverage rate and write results
#         for onebin in sorted_bins:
#             bin_contigs = bin_contig_dir[onebin]
#             n = sum(1 for a in bin_contigs if a in existing_contig_list)
#             cov_rate = n / len(bin_contigs)

#             if cov_rate < 0.2:
#                 existing_contig_list.extend(bin_contigs)
#                 best_bins.append(onebin)

#                 for eachcontig in bin_contigs:
#                     out_contigs2binFile.write(f"{eachcontig}\t{onebin}\n")

#                 # Write bin quality info to output file
#                 out_f.write("\t".join(quality_info_dir[onebin])  + "\n")

#                 # Write contig sequences for each bin to separate FASTA files
#                 with open(os.path.join(output_dir, f"{onebin}.fa"), 'w') as bin_fasta_f:
#                     for contig in bin_contigs:
#                         bin_fasta_f.write(f">{contig}\n{contig_seq_dir[contig]}\n")
######
def main(contig_bin_file, checkm2_file, contig_fasta_file, output_file, output_dir, out_contigs2binF, min_completeness, max_contamination):
    contig_bin_dir = {}         # Mapping: contig ID -> list of candidate bin IDs
    quality_dir = {}            # Mapping: bin ID -> quality score (completeness - 5 * contamination)
    quality_info_dir = {}       # Mapping: bin ID -> original line fields from checkm2
    bin_contig_dir = {}         # Mapping: bin ID -> list of contig IDs
    contig_seq_dir = {}         # Mapping: contig ID -> sequence

    # Parse the input FASTA file and store contig sequences
    parse_fasta(contig_fasta_file, contig_seq_dir)

    # First, parse checkM2 result and contig-bin mapping files
    with open(contig_bin_file, 'r') as contig_bin_f, open(checkm2_file, 'r') as checkm2_f:
        header = checkm2_f.readline()

        # Read and filter bin quality data based on thresholds
        for line in checkm2_f:
            fields = line.strip().split("\t")
            binid = fields[0]
            completeness = float(fields[1])
            contamination = float(fields[2])

            if completeness > min_completeness and contamination < max_contamination:
                qs = completeness - 5 * contamination
                quality_dir[binid] = qs
                quality_info_dir[binid] = fields

        # Read contig-to-bin assignments
        for line in contig_bin_f:
            contigid, binid = line.strip().split("\t")
            contig_bin_dir.setdefault(contigid, []).append(binid)
            bin_contig_dir.setdefault(binid, []).append(contigid)

    # Determine the best-quality bin for each contig
    highest_qs_binids = get_highest_qs_binid(contig_bin_dir, quality_dir)
    highest_qs_binids_QS_dir = {
        binid: quality_dir[binid]
        for binid in highest_qs_binids.values()
        if binid in quality_dir
    }

    # Sort bins in descending order of quality score
    sorted_bins = sorted(highest_qs_binids_QS_dir, key=highest_qs_binids_QS_dir.get, reverse=True)

    # If no bins pass the thresholds, exit early without creating outputs
    if not quality_dir and not sorted_bins:
        print("No valid bins found. Exiting without writing output.")
        return

    # Create output directory only if valid results exist
    create_output_dir(output_dir)

    # Open output files and write results
    with open(output_file, 'w') as out_f, open(out_contigs2binF, 'w') as out_contigs2binFile:
        out_f.write(header)
        best_bins = []
        existing_contig_list = []

        for onebin in sorted_bins:
            bin_contigs = bin_contig_dir[onebin]
            n = sum(1 for a in bin_contigs if a in existing_contig_list)
            cov_rate = n / len(bin_contigs)

            # Skip bins with high contig overlap
            if cov_rate < 0.2:
                existing_contig_list.extend(bin_contigs)
                best_bins.append(onebin)

                # Write contig-to-bin mapping
                for eachcontig in bin_contigs:
                    out_contigs2binFile.write(f"{eachcontig}\t{onebin}\n")

                # Write quality info to report
                out_f.write("\t".join(quality_info_dir[onebin]) + "\n")

                # Write contig sequences to FASTA file for the selected bin
                with open(os.path.join(output_dir, f"{onebin}.fa"), 'w') as bin_fasta_f:
                    for contig in bin_contigs:
                        bin_fasta_f.write(f">{contig}\n{contig_seq_dir[contig]}\n")


######



if __name__ == "__main__":
    # Define argument parser
    parser = argparse.ArgumentParser(description="  Process contig-bin mappings and quality scores to generate bin FASTA files. Author:lianglifeng V1	2023-09-27")
    parser.add_argument('-t','--contig2bin',required=True, help="Path to the contig-bin mapping file.")
    parser.add_argument('-c','--checkm2_file',required=True, help="Path to the CheckM2 results file.")
    parser.add_argument('-f','--contig',required=True, help="Path to the contig FASTA file.")
    parser.add_argument('-r','--report',required=True, help="Path to the output file for bin quality scores.")
    parser.add_argument('-b','--best_bins',required=True, help="Output directory to store bin FASTA files.")
    parser.add_argument('-s','--best_bins_contigs2bin',required=True, help="Output file of best bin contigs2bin.")
    parser.add_argument('-m','--completeness', required=False, type=float, default=50, help="Minimum Completeness threshold (default: 50)")
    parser.add_argument('-n','--contamination', required=False, type=float, default=10, help="Maximum Contamination threshold (default: 10)")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Run the main function
    main(args.contig2bin, args.checkm2_file, args.contig, args.report, args.best_bins, args.best_bins_contigs2bin,args.completeness, args.contamination )
