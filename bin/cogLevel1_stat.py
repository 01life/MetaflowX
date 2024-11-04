#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys


cogCountDir,cog_SYMBOL,cog_DESCRIPTION,cog_CATEGORIE,cogDir={},[],[],[],{}
with open(sys.argv[1],'r') as cogF:

    for a in cogF:
        al = a.strip().split('\t')
        for b in al[0]:
            if b not in cogCountDir:
                cogCountDir[b] = 1
            else:
                cogCount = cogCountDir[b] + 1
                cogCountDir[b] = cogCount

cogCountDF = pd.DataFrame.from_dict(cogCountDir, orient='index', columns=['geneNumber'])
cogLevel = pd.read_csv(sys.argv[2],sep="\t",low_memory=False,index_col=0)
cogLevelCount = pd.concat([cogLevel, cogCountDF], axis=1, join="inner")
cogLevelCount_sorted = cogLevelCount.sort_values(by=["Cog_CATEGORIE", "geneNumber"], ascending=[False, True])
cogLevelCount_sorted['label'] =  cogLevelCount_sorted.index + '|' + cogLevelCount_sorted['Cog_DESCRIPTION']
cogLevelCount_sorted['loggeneNumber'] =  np.log2(cogLevelCount_sorted["geneNumber"])
cogLevelCount_sorted.to_csv(sys.argv[3],sep="\t",na_rep="0")
