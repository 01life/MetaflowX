#!/usr/bin/env python

import os
import argparse
import shutil

def get_fa_files(directory):
    fa_files = {}
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.fa'):
                binID = file.split("_reassembly_contigs")[0]
                fa_files[binID] = os.path.abspath(os.path.join(root, file))
    return fa_files

def filter_and_copy(input_file, fa_dir, hq_dir, lq_dir, lq_bin_id_F, completeness=90, contamination=5):
    fa_files = get_fa_files(fa_dir)
    failed2imporve_bin_list = []
    
    with open(input_file, 'r') as checkmF, open(lq_bin_id_F,'w') as lq_bin_id_File, open('ReAss_Unimproved_Bins_ID.txt','w') as unimproveFile,open('ReAss_improved_info.txt','w') as improveFile :
        header = checkmF.readline()  
        
        for line in checkmF:
            columns = line.strip().split("\t")
            
            binid = columns[0].split("_reassembly_contigs")[0]
            if binid in fa_files:

                # split HQ  or LQ bin fa
                if float(columns[1]) >= completeness and float(columns[2]) <= contamination:
                    dest_path = os.path.join(hq_dir, f"ReAss_HQ_{binid}.fa")
                    improveFile.write(f"{binid}\t{columns[1]}\t{columns[2]}\tsuccessful\tReAss\n")

                    shutil.copy(fa_files[binid], dest_path)
                elif float(columns[1]) <= 50:
                    unimproveFile.write(f"{binid}\n")
                else:
                    lq_bin_id_File.write(str(columns[0])+'\n')
                    # dest_path = os.path.join(lq_dir, fa_filename)
                    shutil.copy(fa_files[binid], lq_dir)
            else:
                print(f'{binid} not in input bins floder')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and copy .fa files based on completeness and contamination thresholds.")
    parser.add_argument('--input_file', help="Input file path")
    parser.add_argument('--fa_dir', help="Directory containing .fa files")
    parser.add_argument('--hq_dir', help="Directory to copy high-quality .fa files")
    parser.add_argument('--lq_bin_id_file', help="LQ bins id list file")
    parser.add_argument('--completeness', type=float, default=90, help="Completeness threshold (default: 90)")
    parser.add_argument('--contamination', type=float, default=5, help="Contamination threshold (default: 5)")
    
    args = parser.parse_args()

    hq_dir = args.hq_dir
    lq_dir = 'LQ_RAB'
    os.makedirs(hq_dir, exist_ok=True)
    os.makedirs(lq_dir, exist_ok=True)

    
    filter_and_copy(args.input_file, args.fa_dir, hq_dir, lq_dir, args.lq_bin_id_file, args.completeness, args.contamination, )
