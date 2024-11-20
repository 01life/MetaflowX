#!/usr/bin/env python

import argparse
import shutil
import os
import glob
from pathlib import Path


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main(deepurify,rename,metainfo,binID):
    deepurify_files = {}
    outPath = Path(rename)
    makedir(outPath)

    for root, _, _ in os.walk(Path(deepurify)):
        files = glob.glob(os.path.join(root, "*.fasta"))

        for file in files:
            deepurify_files[os.path.basename(file)] = file

    binQS = {}
    binQS_DIR = {}

    with open(metainfo, 'r') as metaF:
        for a in metaF:
            al = a.strip().split("\t")
            qs = float(al[1]) - 5 * float(al[2])
            binQS[al[0]] = qs
            binQS_DIR[al[0]] = al[:-1]

    with open(os.path.join(outPath, "deepurify_rename.txt"), 'w') as outFile, open(os.path.join(outPath, "deepurify_rename.QS.txt"), 'w') as newQSF:
        max_bin_name = max(binQS, key=binQS.get)
        newName = f'Deepurify_{binID}.fa'
        shutil.copy(deepurify_files[max_bin_name], os.path.join(outPath, newName))
        outFile.write(f"{newName}\t{max_bin_name}\n")
        qs_txt = "\t".join(binQS_DIR[max_bin_name])
        newQSF.write(f"{binID}\t{qs_txt}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some bin files.")
    parser.add_argument('--deepurify_dir', type=str, help='Path to the deepurify directory')
    parser.add_argument('--rename_dir', type=str, help='Path to the rename directory')
    parser.add_argument('--metaInfo', type=str, help='Path to the MetaInfo.tsv file')
    parser.add_argument('--binid', type=str, help='bin id')
    args = parser.parse_args()
    main(args.deepurify_dir,args.rename_dir,args.metaInfo,args.binid)
