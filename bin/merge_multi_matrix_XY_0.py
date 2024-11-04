#!/usr/bin/env python

import sys
import pandas as pd

uniq_index_ID=[]
profileList=[]
with open(sys.argv[1],'r') as profileF:
    n=0
    for a in profileF:
        al=a.strip().split('\t')
        profile=pd.read_csv(al[1],sep="\t",index_col=0,low_memory=False)
        profile.columns = profile.columns.astype(str)
        profile.index = profile.index.astype(str)
        newColumnsList = [(str(al[0]))]
        profile.columns = newColumnsList
        profileList.append(profile)

result = pd.concat(profileList, axis=1).fillna('0')
result.to_csv(sys.argv[2],sep="\t",na_rep=0)

