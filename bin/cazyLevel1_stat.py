#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

   
   
cazyDir = {}
#cazyCategorieDir = {"GH":"GHs|Glycoside Hydrolases","GT":"GTs|GlycosylTransferases","PL":"PLs|Polysaccharide Lyase","CE":"CEs|Carbohydrate Esterases","AA":"AAs|Auxiliary Activities","CB":"CBMs|Carbohydrate-Binding Modules"}
cazyCategorieDir={}
with open(sys.argv[2],'r') as cataF:
    for b in cataF:
        bl = b.strip().split('\t')
        cazyCategorieDir[bl[0]] = bl[1]


n = 0
with open(sys.argv[1],'r') as cazyF:
    for z in cazyF:
        if n > 0:
            zl = z.strip().split('\t')
            cazyClass = zl[0][:2]
            if cazyCategorieDir[cazyClass] not in cazyDir:
                cazyDir[cazyCategorieDir[cazyClass]] = 1
            else:
                cazyCount = cazyDir[cazyCategorieDir[cazyClass]] + 1 
                cazyDir[cazyCategorieDir[cazyClass]] = cazyCount
        n +=1

# Convert the dictionary to a DataFrame.
cazyCountDF = pd.DataFrame.from_dict(cazyDir, orient='index', columns=['geneNumber'])
cazyCountDF['log2geneNumber'] =  np.log2(cazyCountDF["geneNumber"])
cazyCountDF_sorted = cazyCountDF.sort_values(by=[ "geneNumber"], ascending=[ True])
cazyCountDF_sorted.to_csv(sys.argv[3],sep="\t",na_rep="0")