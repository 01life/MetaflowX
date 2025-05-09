#!/usr/bin/env python

import pandas as pd
import argparse as ap

# MERGEBINPROFILE
def ParsReceiver():
    p = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                          description='''
-------------------
	Funtion:
		Get every contig's GC and avg_depth.
	Version:
		V1	2024-01-10
	Author:
		Wenjing Shi
------------------''')
    ars = p.add_argument
    ars('-d', dest="depthlist", action="store", help='alldepth.list', required=True)
    # # CONTIGSTAT
    ars('-i', dest="contig_info", action="store", help='all_contig_info.txt', required=True)
    # RENAMEBIN
    ars('-b', dest="contig2bin", action="store", help='bin.fa.list', required=True)
    ars('-o', dest="gc_depth", action="store", help='MetaFlowX_contigs_gc_depth.xls', required=True)
    return vars(p.parse_args())


pars = ParsReceiver()

a = pars['depthlist']
b = pars['contig_info']
c = pars['contig2bin']
d = pars['gc_depth']

df1 = pd.read_table(a, header=None)
i = 0
final_df = pd.DataFrame()
while i < len(df1):
    sample = df1.iloc[i, 0]
    dep_path = df1.iloc[i, 1]
    dep_df = pd.read_table(dep_path, usecols=['totalAvgDepth'])
    dep_df.rename(columns={'totalAvgDepth': sample}, inplace=True)
    final_df = pd.concat([final_df, dep_df], axis=1)
    i += 1
final_df['avg_depth'] = final_df.mean(axis=1)
bins_order = pd.read_table(dep_path, usecols=['contigName'])
new_df1 = pd.concat([bins_order, final_df], axis=1)
# 20241213 由于contigName内容增加，修改
new_df1['contigName'] = new_df1['contigName'].apply(lambda x: "|".join(x.split("|")[:2]))
##

df2 = pd.read_table(b)
new_df2 = pd.merge(new_df1, df2, left_on='contigName', right_on='ContigID', how='left')
with open(c, 'r') as inputF:
    binlist = []
    contiglist = []
    for a in inputF:
        al = a.strip().split('\t')
        with open(al[1], 'r') as binF:
            for line in binF:
                if line.startswith('>'):
                    binlist.append(al[0])
                    contiglist.append(line[1:-1])
df3 = pd.DataFrame({'BinID': binlist, 'contigName': contiglist})
# 20241213 由于contigName内容增加，修改
df3['contigName'] = df3['contigName'].apply(lambda x: "|".join(x.split("|")[:2]))
##
final_df = pd.merge(df3, new_df2, on='contigName', how='left')
final_df = final_df[['BinID', 'contigName', 'avg_depth', 'GC']]
final_df.to_csv(d, sep='\t', index=False)
