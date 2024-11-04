#!/usr/bin/env python

import argparse
import shutil
from pathlib import Path
from typing import Dict
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
        binid = file.stem.split("Deepurify_")[1]
        fa_files[binid] = str(file.resolve())
    return fa_files


def file2IDir(oneFile):
    tmpDir = {}
    with open(oneFile,'r') as oneF:
        for a in oneF:
            al = a.strip().split("\t")
            tmpDir[al[0]] = al[1:]
    return(tmpDir)



def main(mergeQSFile: str,deepurify_infoFile:str, deepurify_result: str,   outfile: str, improveFile:str):
    """
    Main function to process the input file and copy the best quality bin.
    
    :param input_file: Path to the input TSV file
    :param deepurify_result: Path to the Deepurify result directory
    :param bin_id: Bin ID to be used in output file names
    :param org_Completeness: Original completeness threshold
    :param org_Contamination: Original contamination threshold
    """
    # Input validation
    if not Path(mergeQSFile).is_file():
        raise FileNotFoundError(f"Input file not found: {mergeQSFile}")
    if not Path(deepurify_infoFile).is_file():
        raise FileNotFoundError(f"Input file not found: {mergeQSFile}")

    if not Path(deepurify_result).is_dir():
        raise NotADirectoryError(f"DAS Tools result directory not found: {deepurify_result}")

    org_rebin_qs_dir = file2IDir(mergeQSFile)
    deepurify_qs_dir = file2IDir(deepurify_infoFile)

    deepurify_fa_dir = get_fa_files(deepurify_result)

    failed_improve_bin = []
    improve_info = []
    for org_binid in deepurify_qs_dir:
        # org_binid = onebin.strip().lstrip("Deepurify_").rstrip("\.fa")
        deepurify_com = float(deepurify_qs_dir[org_binid][1])
        deepurify_cov = float(deepurify_qs_dir[org_binid][2])
        deepurify_QS = deepurify_com - 5*deepurify_cov
        if org_binid in org_rebin_qs_dir:
            org_com,org_cov,rebin_com,rebin_cov = org_rebin_qs_dir[org_binid][1:]
            org_QS = float(org_com) - 5* float(org_cov)

        else:
            print(f"{org_binid} not in deepurify_qs_dir ")

        pick_bin = ''
        

        if deepurify_com > float(org_com) and deepurify_cov < float(org_cov):
            pick_bin = deepurify_fa_dir[org_binid]
            output_file = f'ReAss_ReBin_ReRefine_HQ_{org_binid}.fa'
            shutil.copy(pick_bin, output_file)
            improve_info.append( f"{org_binid}\t{deepurify_com}\t{deepurify_cov}\tsuccessful\tRefine" )

        elif deepurify_QS > org_QS:
            pick_bin = deepurify_fa_dir[org_binid]
            output_file = f'ReAss_ReBin_ReRefine_Max_QS_{org_binid}.fa'
            shutil.copy(pick_bin, output_file)
            improve_info.append( f"{org_binid}\t{deepurify_com}\t{deepurify_cov}\tsuccessful\tRefine")
        else:
            failed_improve_bin.append(org_binid)

    with open(outfile,'w') as failedFile:
        if failed_improve_bin != [] :
            failedFile.write('\n'.join(failed_improve_bin)+'\n')

    with open(improveFile,'w') as improveF:   
        if improve_info != []:
            improveF.write('\n'.join(improve_info))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Based on the re-binning checkm result, pick the best quality rebinning bin.')
    parser.add_argument('--mergeQS', type=str, help='meger original and re-binning checkm2 result')
    parser.add_argument('--deepurify_info', type=str, help='deepurify_rename.QS.txt')
    parser.add_argument('--deepurify_fa', type=str, help='path of Deepurify_Result')
    parser.add_argument('--failed', type=str, help='Failed Quality Improvement Reassembly Bin txt')
    parser.add_argument('--improve', type=str, help='improve Bin info txt')
    args = parser.parse_args()

    for arg_name, arg_value in vars(args).items():
        if arg_value is None:
            parser.error(f"Missing required argument: --{arg_name}")

    # try:
    main(args.mergeQS, args.deepurify_info, args.deepurify_fa, args.failed, args.improve)
    # except ValueError as e:
    #     logging.error(f"Error: {e}")
    # except Exception as e:
    #     logging.error(f"An unexpected error occurred: {e}")