#!/usr/bin/env python

import sys,Bio
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import shutil
from pathlib import Path
import os
import glob

def makedir(path):
	if not os.path.exists(path):
		os.makedirs(path)

def most_common_element(binID_list):
    # Use Counter to count the occurrences of each element in the list.
    counts = Counter(binID_list)
    # Use the most_common method to find the most frequently occurring elements and their counts.
    most_common_item, most_common_count = counts.most_common(1)[0]
    return most_common_item


def get_contig_bin(binID,binfa, cDir):
    for record in SeqIO.parse(binfa, 'fasta'):
        cDir[record.id] = binID

deepurify_files,contig_bin,newNamePath = {},{},{}


# Step 1: Get .fa / fasta files from the filter bin directory.
for root, _, _ in os.walk(Path(sys.argv[1])):
    fa_files = glob.glob(os.path.join(root, "*.fa"))
    fasta_files = glob.glob(os.path.join(root, "*.fasta"))
    files = fa_files + fasta_files

    for file in files:
        file_name = os.path.splitext(os.path.basename(file))[0]            
        get_contig_bin(file_name,file,contig_bin)
            
most_bin_dir={}

outPath = Path(sys.argv[3])
makedir(outPath)
# Step 2: Get .fa / fasta files from the deepurify directory.
for root, _, _ in os.walk(Path(sys.argv[2])):
    fa_files = glob.glob(os.path.join(root, "*.fa"))
    fasta_files = glob.glob(os.path.join(root, "*.fasta"))
    files = fa_files + fasta_files
    for file in files:
        file_name = os.path.splitext(os.path.basename(file))[0].lstrip("Deepurify_").rstrip(".fasta").rstrip(".fa")
        deepurify_files[file_name] = file
        binID_list= []
        for record in SeqIO.parse(file, 'fasta'):
            binID_list.append(contig_bin[record.id])
        most_bin = most_common_element(binID_list)
        most_bin_dir.setdefault(most_bin,[]).append(file_name)

binQS = {}
binQS_DIR = {}
with open(sys.argv[4],'r')  as metaF:
    for a in metaF:
        al =a.strip().split("\t")
        qs = float(al[1]) - 5*float(al[2])
        binID_s = al[0].lstrip("Deepurify_").rstrip(".fasta").rstrip(".fa")
        binQS[binID_s]= qs
        binQS_DIR[binID_s]= al[:-1]

with open(os.path.join(outPath,"deepurify_rename.txt"),'w') as outFile,open(os.path.join(outPath,"deepurify_rename.QS.txt"),'w') as newQSF:

    for b in most_bin_dir:
        newName = f'Deepurify_{b}.fa'
        if len(most_bin_dir[b]) == 1:
            shutil.copy(deepurify_files[most_bin_dir[b][0]], os.path.join(outPath, newName))
            outFile.write(f"{most_bin_dir[b][0]}\t{newName}\n")
            qs_txt = "\t".join(binQS_DIR[most_bin_dir[b][0]])
            newQSF.write(f"{newName}\t{qs_txt}\n")
            
        else:
            # Get the corresponding qs values for each bin in the current bin list.
            bin_qs_values = {bin_name: binQS.get(bin_name, 0) for bin_name in most_bin_dir[b]}
            # Find the bin name with the maximum qs value.
            max_bin_name = max(bin_qs_values, key=bin_qs_values.get)
            shutil.copy(deepurify_files[max_bin_name], os.path.join(outPath, newName))
            outFile.write(f"{max_bin_name}\t{newName}\n")
            qs_txt = "\t".join(binQS_DIR[max_bin_name])
            newQSF.write(f"{newName}\t{qs_txt}\n")



