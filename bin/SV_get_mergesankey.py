#!/usr/bin/env python

import pandas as pd 
import sys

input1 = pd.read_table(sys.argv[1])
input2 = pd.read_table(sys.argv[2])
input1.columns = input1.columns.str.replace('_bin_abundance', '')

input1['sum'] = input1.iloc[:, 1:].sum(axis=1)
input1 = input1.sort_values(by='sum', ascending=False)

input2 = input2.drop(input2.columns[1], axis=1)
df = pd.merge(input2, input1, how='right')
df.to_csv(sys.argv[3], sep='\t', index=False)