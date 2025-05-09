#!/usr/bin/env python

import pandas as pd
import argparse as ap
import os


# MERGEBINPROFILE
def ParsReceiver():
    p = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                          description='''
-------------------
	Funtion:
		Merge bins species annotations and abundance files.
	Version:
		V1	2024-01-09
	Author:
		Wenjing Shi
------------------''')
    ars = p.add_argument
    ars('-i', dest="anno", action="store", help='MetaFlowX_bins_taxonomic.xls', required=True)
    ars('-a', dest="abundancefile", help='MetaFlowX_bins_TPM.xls', action='store', required=True)
    # output
    ars('-o', dest="outpath", action="store", help='output path', required=True)
    ars('-p', dest="prefix", action="store", help='{prefix}_bins_{level}.xls', required=True)
    ars('-l', dest="level", action="store", help='Taxon[species, genus, family, order, class, phylum]', required=True)
    return vars(p.parse_args())


pars = ParsReceiver()
df1 = pd.read_table(pars['anno'])
df2 = pd.read_table(pars['abundancefile'])
df2.columns = [col.replace("_bin_abundance", "") for col in df2.columns]

def merge(taxon):
    global df1, df2
    taxon_df = df1[['BinID', taxon]].copy()
    df = pd.merge(taxon_df, df2, on=['BinID'])
    df = df.drop('BinID', axis=1)
    new_df = df.groupby(taxon).sum()
    return new_df


i = pars['level'].capitalize()
if i in df1.columns[2:8]:
    path = os.path.join(os.path.normpath(pars['outpath']), pars['prefix'] + '_bins_' + i.lower() + '.xls')
    merge(i).to_csv(path, sep='\t')


# if len(sys.argv) == 1:
# 	sys.argv.append('-h')
# 	ParsReceiver()
# 	sys.exit()
