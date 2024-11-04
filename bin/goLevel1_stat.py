#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

oneCategorieDir, oneDir, oneCountDir={},{},{}
with open(sys.argv[1],'r') as oneProfileF, open(sys.argv[2],'r') as oneDBcategoryF:
    for o in oneDBcategoryF:
        ol = o.strip().split('\t')
        oneCategorieDir[ol[0]] = ol[1]
        if ol[1] not in oneCountDir:
            oneCountDir[ol[1]] = 0
    for i in oneProfileF:
        il = i.strip().split('\t')
        if il[0] in oneCategorieDir:
            oneCount =  oneCountDir[oneCategorieDir[il[0]]] + 1
            oneCountDir[oneCategorieDir[il[0]]] = oneCount
        else:
            if 'Others' not in oneCountDir:
                oneCountDir['Others'] = 1
            else:
                oneCount =  oneCountDir['Others'] + 1
                oneCountDir['Others'] = oneCount
oneCountDF = pd.DataFrame.from_dict(oneCountDir, orient='index', columns=['geneNumber'])
oneCountDF['log2geneNumber'] =  np.log2(oneCountDF["geneNumber"])
oneCountDF_sorted = oneCountDF.sort_values(by=[ "geneNumber"], ascending=[ True])

oneCountDF_sorted.to_csv(sys.argv[3],sep="\t",na_rep="0")
