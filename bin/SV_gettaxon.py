#!/usr/bin/env python

import sys
import pandas as pd

input = sys.argv[1]
output = sys.argv[2]
# ADDNCBI
df = pd.read_table(input, usecols=['user_genome', 'classification'])
new_df = df['classification'].str.split(';', expand=True)
new_df = new_df.apply(lambda col: col.str[3:])
df = pd.concat([df, new_df], axis=1)
df = df.drop('classification', axis=1)
df['index'] = df['user_genome'].str[4:].astype(int)
df = df.sort_values(by='index').reset_index(drop=True)
# unclassified species
i = 0
while i < len(df):
    if len(df.loc[i, 6]) == 0:
        df.loc[i, 6] = df.loc[i, 5] + ' uc_sp'
    i += 1
df = df.drop('index', axis=1)
df.columns = ['BinID', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
df['Species'] = df['Species'].str.replace(' ', '_')
df.to_csv(output, index=False, sep='\t')
