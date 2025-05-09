#!/usr/bin/env python

import sys
import pandas as pd

input = sys.argv[1]
output = sys.argv[2]

df = pd.read_table(input, sep='\t', index_col=0)
# The format is the same as TPM, RPKM, and AvgDepth.
# df.columns = df.columns.to_series().replace(r'\.sorted.*$', '_bin_abundance', regex=True)
# df = df.iloc[1:]

df = df[df.index != 'unmapped']
df = df.div(df.sum(axis=0), axis=1) * 100
df = df.rename_axis('BinID')
df.to_csv(output, index=True, sep='\t')
