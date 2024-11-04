#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def create_output_dir(output_dir):
    """
    Create output directory if it doesn't exist.
    """
    os.makedirs(output_dir, exist_ok=True)

def rename_fasta(old_fa, new_fa, binID):
    """
    rename contig id of one bin.
    """
    with open(new_fa,'w') as new_faste:
        for record in SeqIO.parse(old_fa, "fasta"):
            newID = f"{record.id}|{binID}"
            new_faste.write(f">{newID}\n{record.seq}\n")


def rename_files_and_log(input_dir, output_file_prefix, output_dir):
    n = 1
    rename_map_file = f'{output_file_prefix}_HQ_unique_bins_rename_map.xls'

    with open(rename_map_file, 'w') as map_file:
        for file in os.listdir(input_dir):
            if file.endswith('.fa') or file.endswith('.fasta') or file.endswith('.fna'):  
                file_path = os.path.join(input_dir, file)
                if os.path.isfile(file_path):
                    newname = f'bin.{n}.fa'
                    newname_path = os.path.join(output_dir, newname)
                    rename_fasta(file_path, newname_path,f'bin.{n}')
                    map_file.write(f'{file} {newname}\n')
                    n += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename files and log the changes.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing the files to rename.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='rename bin floder.')
    parser.add_argument('-p', '--pipeline_prefix', type=str, required=True, help='Prefix for the output log file.')
    
    args = parser.parse_args()

    create_output_dir(args.output_dir)

    rename_files_and_log(args.input_dir, args.pipeline_prefix, args.output_dir)
