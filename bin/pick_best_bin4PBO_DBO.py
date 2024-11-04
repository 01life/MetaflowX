#!/usr/bin/env python

import os
import shutil
import argparse
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Pick the best bins result from PBO and DBO methods. Author: Lianglifeng V1 2023-09-27")
    parser.add_argument('-pb', '--pbo_bins', required=True, help="PBO bins folder")
    parser.add_argument('-ps', '--pbo_quality', required=True, help="PBO checkm2 report")
    parser.add_argument('-bd', '--dbo_bins', required=True, help="DBO bins folder")
    parser.add_argument('-bs', '--dbo_quality', required=True, help="DBO checkm2 report")
    parser.add_argument('-s', '--sample', required=True, help="Sample ID")
    parser.add_argument('-a','--completeness', required=False, type=float, default=50, help="Minimum Completeness threshold (default: 50)")
    parser.add_argument('-b','--contamination', required=False, type=float, default=10, help="Maximum Contamination threshold (default: 10)")
    parser.add_argument('-r','--similarity_ratio', required=False, type=float, default=0.8, help="similarity bins threshold (default: 0.8)")

    return parser.parse_args()

def create_output_dir(output_dir):
    """
    Create output directory if it doesn't exist.
    """
    os.makedirs(output_dir, exist_ok=True)

def parse_single_fasta(fasta_path):
    """
    Parse a single fasta file and return the file name and list of contigs.
    """
    contig_list = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_list.append(record.id)
    return os.path.basename(fasta_path), contig_list

def parse_fasta_directory(fasta_dir):
    """
    Parse all fasta files in the directory using parallel processing.
    """
    bin_contig_dir = {}
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fa") or f.endswith(".fasta")]

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(parse_single_fasta, fasta_file): fasta_file for fasta_file in fasta_files}

        for future in as_completed(futures):
            binID, contigs = future.result()
            bin_contig_dir[binID] = contigs

    return bin_contig_dir

def parse_quality_file(quality_file, min_com, max_con):
    """
    Parse the quality file and return a dictionary of bin quality scores.
    """
    bin_qs_dir = {}
    bin_info_dir = {}
    with open(quality_file, 'r') as f:
        header = f.readline()  # Skip header line
        for line in f:
            data = line.strip().split("\t")
            binid = f"{data[0]}.fa"
            completeness = float(data[1])
            contamination = float(data[2])
            if completeness > min_com  and contamination < max_con:
                qs = completeness - 5 * contamination
                bin_qs_dir[binid] = qs  # data[0] is the bin ID
                bin_info_dir[binid] = data
    return bin_qs_dir,bin_info_dir,header


def pick_best_bin(pbo_qs_dir, dbo_qs_dir, pbo_contig_dir, dbo_contig_dir, pbo_bins_dir, dbo_bins_dir,outpath,pbo_bin_info_dir,dbo_bin_info_dir,qs_report,qs_header,cov_result,similarity_ratio):

    """
    Pick the best bins based on quality scores and contig overlap between PBO and DBO methods.
    """
    hq_bins = []
    similar_pbo_bins = []
    # cov_dir = {}


    cov_result_out = open(cov_result,'w')
    cov_result_out.write("binA\tbinB\tratio\n")

    for dbo_bin, dbo_qs in dbo_qs_dir.items():
        dbo_contigs = dbo_contig_dir[dbo_bin]
        tmp_bin_qs = {dbo_bin: dbo_qs}

        for pbo_bin, pbo_qs in pbo_qs_dir.items():
            pbo_contigs = pbo_contig_dir[pbo_bin]

            # Calculate coverage rate
            intersection = set(pbo_contigs).intersection(set(dbo_contigs))
            cov_rate = len(intersection) / len(dbo_contigs)
            
            if cov_rate > similarity_ratio:
                tmp_bin_qs[pbo_bin] = pbo_qs
                similar_pbo_bins.append(pbo_bin)
                
            # if cov_rate > 0:
            #     cov_dir.setdefault(dbo_bin,{})[pbo_bin] = cov_rate
            cov_result_out.write(f"{dbo_bin}\t{pbo_bin}\t{cov_rate}\n")

        # Pick the bin with the highest quality score
        highest_qs_binid = max(tmp_bin_qs, key=tmp_bin_qs.get)
        hq_bins.append(highest_qs_binid)


    # Add remaining PBO bins that were not similar
    alone_pbo = set(pbo_qs_dir.keys()) - set(similar_pbo_bins)

    for a_pbo_bin in alone_pbo:
        tmp_pbo_qs = {a_pbo_bin: pbo_qs_dir[a_pbo_bin]}

        odd_alone_pbo_list = [item for item in alone_pbo if item != a_pbo_bin]
        a_pbo_bin_contigs = pbo_contig_dir[a_pbo_bin]


        for odd_pbo in odd_alone_pbo_list:
            odd_contigs = pbo_contig_dir[odd_pbo]
            pbo_intersection = set(a_pbo_bin_contigs).intersection(set(odd_contigs))
            a_pbo_cov_rate = len(pbo_intersection) / len(a_pbo_bin_contigs)
            cov_result_out.write(f"{a_pbo_bin}\t{odd_pbo}\t{a_pbo_cov_rate}\n")

            if a_pbo_cov_rate > similarity_ratio:
                tmp_pbo_qs[odd_pbo] = pbo_qs_dir[odd_pbo]

        highest_qs_pbo = a_pbo_bin
        if len(tmp_pbo_qs.keys()) > 1 :
            highest_qs_pbo = max(tmp_pbo_qs, key=tmp_pbo_qs.get)

        if highest_qs_pbo not in hq_bins:
            hq_bins.append(a_pbo_bin)


    # Copy corresponding fasta files to Best_Bin folder
    with open(qs_report,'w') as qs_reportFile:
        qs_reportFile.write(str(qs_header.strip())+"\n")
        for bin_id in set(hq_bins):
            pbo_bin_file = os.path.join(pbo_bins_dir, bin_id )
            dbo_bin_file = os.path.join(dbo_bins_dir, bin_id )

            # Copy if file exists in PBO bins
            if os.path.exists(pbo_bin_file):
                shutil.copy(pbo_bin_file, outpath)
            # Copy if file exists in DBO bins
            elif os.path.exists(dbo_bin_file):
                shutil.copy(dbo_bin_file, outpath)
            
            out_bin_id = bin_id.replace(".fa","")

            if bin_id in pbo_bin_info_dir:
                qs_reportFile.write(str(out_bin_id)+"\t"+"\t".join(pbo_bin_info_dir[bin_id][1:])+'\n')
            
            elif bin_id in dbo_bin_info_dir:
                qs_reportFile.write(str(out_bin_id)+"\t"+"\t".join(dbo_bin_info_dir[bin_id][1:])+'\n')

    best_bin_contig_dir = parse_fasta_directory(outpath)
    best_contig_bin_dir = {}

    for bb, bb_contigs in best_bin_contig_dir.items():
        for a_contig in bb_contigs:
            best_contig_bin_dir.setdefault(a_contig,[]).append(bb)

    cov_result_out.close()


def main():
    args = parse_args()


    outpath = f"{args.sample}_Best_Bin"
    create_output_dir(outpath)

    # Parse PBO and DBO bin fasta directories
    pbo_bin_contig_dir = parse_fasta_directory(args.pbo_bins)
    dbo_bin_contig_dir = parse_fasta_directory(args.dbo_bins)

    # Parse PBO and DBO quality files
    min_com = args.completeness 
    max_con = args.contamination
    pbo_bin_qs_dir,pbo_bin_info_dir,pbo_qs_header = parse_quality_file(args.pbo_quality,min_com, max_con)
    dbo_bin_qs_dir,dbo_bin_info_dir,dbo_qs_header = parse_quality_file(args.dbo_quality,min_com, max_con )

    # Pick the best bins
    qs_report = f"{args.sample}_Best_Bin_quality_report.tsv"
    cov_result = f"{args.sample}_bins_contig_similarity_ratio.csv"
    pick_best_bin(pbo_bin_qs_dir, dbo_bin_qs_dir, pbo_bin_contig_dir, dbo_bin_contig_dir,args.pbo_bins,args.dbo_bins,outpath,pbo_bin_info_dir, dbo_bin_info_dir,qs_report,dbo_qs_header,cov_result, args.similarity_ratio)

if __name__ == "__main__":
    main()
