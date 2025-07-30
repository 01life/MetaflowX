#!/usr/bin/env python

import sys
import os
from pathlib import Path


deepurify_path = Path(sys.argv[1]).resolve()
deepurify_files = {}

for root, _, files in os.walk(deepurify_path):
    for file in files:
        if file.endswith('.fa'):
            file_path = os.path.join(root, file)
            file_name = os.path.splitext(file)[0].lstrip("Deepurify_")
            deepurify_files[file_name] = file_path


with open(sys.argv[2],'r')  as orgmapF, open(sys.argv[3],'w') as outF:
    header = orgmapF.readline()
    outF.write(header)
    for a in orgmapF:
        binID,binGenome,refGenomes,highestAbundanceFq1,fq1 = a.strip().split("\t")
        if binID in deepurify_files:
            binGenome = deepurify_files[binID]
            outF.write(f"{binID}\t{binGenome}\t{refGenomes}\t{highestAbundanceFq1}\t{fq1}\n")