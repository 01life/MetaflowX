#!/usr/bin/env python

import pandas as pd 
import argparse as ap


def ParsReceiver():
    p = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                          description='''
-------------------
	Funtion:
		Get tree merge file.
	Version:
		V1	2024-03-19
	Author:
		Wenjing Shi
------------------''')
    ars = p.add_argument
    ars('-i', dest="input", action="store", help='tree file [.tree/.nwk]', required=True)
    ars('-t', dest="taxonomic", help='MetaFlowX_bins_taxonomic.xls', action='store', required=True)
    ars('-a', dest="info", help='MetaFlowX_bins_info.xls', action='store', required=True)
    ars('-o', dest="output", help=' ', action='store', required=True)
    return vars(p.parse_args())


pars = ParsReceiver()

df1 = pd.read_table(pars['taxonomic'])
df2 = pd.read_table(pars['info'])
out = pars['output']

df1 = df1.iloc[:, :3]
df = pd.merge(df1, df2, how='left')

#
with open(pars['input'], "r") as tree_file:
    tree_data = tree_file.read()
    tree_file.close()
df.to_csv(out, sep='\t', index=False)
with open(out, "a") as output_file:
    output_file.write('#\n'+tree_data)