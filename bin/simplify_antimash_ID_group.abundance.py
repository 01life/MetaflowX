#!/usr/bin/env python

import pandas as pd
import re,sys


df = pd.read_csv(sys.argv[1],index_col=0)

colIndexList = ['RAW','RPKM','TPM','coreRAW','coreRPKM','coreTPM','cov','corecov']
columnIndexDir = {}

for a in df.columns.tolist():
    atype=a.strip().split('.')[-1]
    columnIndexDir[atype] = a


df['newID'] = df.index.map(lambda x: re.findall(r'Entryname=(.*?)--', str(x))[0])

sampleName = sys.argv[2].strip()


for b in colIndexList:
    if b in columnIndexDir:
        grouped_raw = df.groupby('newID')[columnIndexDir[b]].sum()
        grouped_raw = grouped_raw.rename(str(sampleName))
        grouped_raw.to_csv(sys.argv[3]+'_'+str(b)+'.xls',sep="\t",na_rep="NA",header=True)

