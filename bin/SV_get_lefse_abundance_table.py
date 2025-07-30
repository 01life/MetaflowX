#!/usr/bin/env python

import argparse as ap
import pandas as pd


def ParsReceiver():
    p = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                          description='''
-------------------
	Funtion:
		Get abundance_table for lefse.
	Version:
		V1	2024-02-19
	Author:
		Wenjing Shi
------------------''')
    ars = p.add_argument
    ars('-i', dest="input", action="store", help='gtdbtk.bac120.summary.tsv', required=True)
    ars('-a', dest="abundance", help='MetaFlowX_CoverM_bins_relative_abundance.xls', action='store', required=True)
    ars('-o', dest="output", help='MetaFlowX_bins_abundance_table.xls', action='store', required=True)
    return vars(p.parse_args())


pars = ParsReceiver()
input = pars['input']
ab = pars['abundance']
output = pars['output']

df1 = pd.read_table(input, usecols={'user_genome', 'classification'})
df2 = pd.read_table(ab)
df2.columns = [col.replace("_bin_abundance", "") for col in df2.columns]
df = pd.merge(df1, df2, left_on='user_genome', right_on = 'BinID')
df = df.drop(df.columns[[0, 2]], axis=1)
df['classification'] = df['classification'].str.replace(';', '|').str.replace('d__', 'k__').str.replace(' ', '_')

final_df = pd.DataFrame()
for i in range(1,8):
    new_df = df.copy()
    new_df['clade_name'] = new_df['classification'].str.split('|').str[:i].str.join('|')
    new_df= new_df.groupby('clade_name').sum(numeric_only=True)
    new_df['clade_name'] = new_df.index
    final_df = pd.concat([final_df, new_df], ignore_index=True)

final_df = final_df[final_df.columns.tolist()[-1:] + final_df.columns.tolist()[:-1]]
final_df.to_csv(output, sep='\t', index=False)