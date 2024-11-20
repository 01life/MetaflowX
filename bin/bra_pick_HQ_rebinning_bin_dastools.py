#!/usr/bin/env python

import pandas as pd
import argparse
import shutil
from pathlib import Path
from typing import Dict, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_fa_files(directory: str) -> Dict[str, str]:
    """
    Recursively find all .fa files in the given directory.
    
    :param directory: Path to the directory to search
    :return: Dictionary with file stems as keys and full paths as values
    """
    fa_files = {}
    for file in Path(directory).rglob('*.fa'):
        fa_files[file.stem] = str(file.resolve())
    return fa_files

def get_org_quality(org_F: str, binID: str) -> Tuple[float, float]:
    """
    Get the completeness and contamination for a given bin ID.
    
    :param org_F: Path to the input file
    :param binID: ID of the bin to search for
    :return: Tuple of (completeness, contamination)
    :raises ValueError: If the bin ID is not found
    """
    with open(org_F, 'r') as org_File:
        header = org_File.readline()  # ignore header
        for a in org_File:
            al = a.strip().split('\t')
            if al[0] == binID:
                return float(al[2]), float(al[3])
    raise ValueError(f"No matching bin ID found: {binID}")

def main(input_file: str, das_tools_result: str, bin_id: str, org_Completeness: float, org_Contamination: float, outfile: str, pick_checkm:str):
    """
    Main function to process the input file and copy the best quality bin.
    
    :param input_file: Path to the input TSV file
    :param das_tools_result: Path to the DAS Tools result directory
    :param bin_id: Bin ID to be used in output file names
    :param org_Completeness: Original completeness threshold
    :param org_Contamination: Original contamination threshold
    """

    improve_info=''
    # Input validation
    if not Path(input_file).is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if not Path(das_tools_result).is_dir():
        raise NotADirectoryError(f"DAS Tools to checkm2 result directory not found: {das_tools_result}")

    df = pd.read_csv(input_file, sep="\t", low_memory=False)
    bin_fa = get_fa_files(das_tools_result)

    filtered_df = df[(df['Completeness'] > org_Completeness) & (df['Contamination'] < org_Contamination)].copy()


    if not filtered_df.empty:
        filtered_df.loc[:, 'QS'] = filtered_df['Completeness'] - 5 * filtered_df['Contamination']
        best_max_qs_id = filtered_df.loc[filtered_df['QS'].idxmax(), 'Name']
        output_file = f'ReAss_ReBin_HQ_{bin_id}.fa'
        improve_info += f"{bin_id}"
    else:
        df.loc[:, 'QS'] = df['Completeness'] - 5 * df['Contamination']
        best_max_qs_id = df.loc[df['QS'].idxmax(), 'Name']
        if df.loc[df['Name'] == best_max_qs_id, 'Completeness'].values[0] <= 50:
            output_file = f'ReAss_ReBin_Unimprove_{bin_id}.fa'
        else:
            output_file = f'ReAss_ReBin_Max_QS_{bin_id}.fa'


    best_completeness = df.loc[df['Name'] == best_max_qs_id, 'Completeness'].values[0]
    best_contamination = df.loc[df['Name'] == best_max_qs_id, 'Contamination'].values[0]

    # Get all information where the Name column equals best_max_qs_id.
    pick_df = df.loc[df['Name'] == best_max_qs_id].copy()
    # Write this information to a file, such as a CSV file.
    pick_df['Name'] = [bin_id]
    pick_df.to_csv(pick_checkm, index=False,sep="\t")


    with open(outfile,'w') as outF:
        outF.write(f'{bin_id}\t{output_file}\t{org_Completeness}\t{org_Contamination}\t{best_completeness}\t{best_contamination}\n')
    
    with open(outfile+'_improve_info.txt','w') as outF2:
        if improve_info != "":
            outF2.write(f'{improve_info}\t{best_completeness}\t{best_contamination}\tsuccessful\tReBin\n')

    try:
        shutil.copy(bin_fa[best_max_qs_id], output_file)
        logging.info(f"File successfully copied to {output_file}")
    except IOError as e:
        logging.error(f"An error occurred while copying the file: {e}")
    except KeyError as e:
        logging.error(f"The specified bin ID was not found: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Based on the re-binning checkm result, pick the best quality rebinning bin.')
    parser.add_argument('--rebin_checkm', type=str, help='Path to the input TSV file')
    parser.add_argument('--das_tools_result', type=str, help='Path to the DAS Tools result directory')
    parser.add_argument('--bin_id', type=str, help='Bin ID to be used in output file names')
    parser.add_argument('--org_checkm', type=str, help='bin_QS_taxonomy_summary.xls')
    parser.add_argument('--mergeQS', type=str, help='meger original and re-binning checkm2 result')
    parser.add_argument('--pick_checkm', type=str, help='the checkm2 result of picking re-binning bin')

    args = parser.parse_args()

    for arg_name, arg_value in vars(args).items():
        if arg_value is None:
            parser.error(f"Missing required argument: --{arg_name}")

    try:
        org_Completeness, org_Contamination = get_org_quality(args.org_checkm, args.bin_id)
        main(args.rebin_checkm, args.das_tools_result, args.bin_id, org_Completeness, org_Contamination,args.mergeQS, args.pick_checkm)
    except ValueError as e:
        logging.error(f"Error: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")