#!/usr/bin/env python

import argparse as ap
import pandas as pd

# ADDNCBI
def ParsReceiver():
    p = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                          description='''
-------------------
	Funtion:
		Get final bins info file.
	Version:
		V1	2024-01-09
	Author:
		Wenjing Shi
------------------''')
    ars = p.add_argument
    ars('-i', dest="info1", action="store", help='MetaFlowX_raw_bin_quality.xls', required=True)
    ars('-r', dest="rename", help='MetaFlowX_All_bins_rename_map.xls', action='store', required=True)
    ars('-b', dest="info2", action="store", help='MetaFlowX_All_bins_info.xls', required=True)
    ars('-t', dest="taxonomic", help='MetaFlowX_bins_taxonomic.xls', action='store', required=True)
    ars('-o', dest="output", help='MetaFlowX_bins_info.xls', action='store', required=True)
    return vars(p.parse_args())


pars = ParsReceiver()
input = pars['info1']
rename = pars['rename']
info = pars['info2']
taxonomic = pars['taxonomic']
output = pars['output']

df1 = pd.read_table(input)
df2 = pd.read_table(rename, header=None, sep=' ')

df2[0] = df2[0].str[:-3]
dict = df2.set_index(0)[1].to_dict()
df1['BinID'] = df1['Name'].map(dict).str[:-3]
newdf1 = df1
df3 = pd.read_table(info)
newdf2 = pd.merge(newdf1, df3, on='BinID', how='left')
newdf2 = newdf2[['BinID', 'contigNumber', 'GenomeSize', 'Completeness', 'QS', 'Contamination', 'Contig_N50', 'GC_Content']]

df4 = pd.read_table(taxonomic)
df4 = df4[['BinID', 'Species']]
final_df = pd.merge(newdf2, df4, on='BinID', how='left')
# Use fillna; otherwise, a warning will occur during to_csv (due to missing values), but this data is needed when plotting.
final_df.fillna('NotBin', inplace=True)
final_df.to_csv(output, sep='\t', index=False)
