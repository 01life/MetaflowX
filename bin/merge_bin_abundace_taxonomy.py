#!/usr/bin/env python

import pandas as pd
import argparse as ap


def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    Function:
        Merge bin information to show in report
    
    Version:
        V1 2024-05-29
    Author:
        Lifeng Liang
-------------------'''
)
    ars = p.add_argument
    ars("-p",dest='profile_file', help='Bin abundance profile file')
    ars("-t",dest='taxonomy_file', help='Bin GTDB taxonomy file')
    ars("-s",dest='show_unmap', choices=['Y', 'N'], help='Flag to include unmapped entries',default='Y')
    ars('-q',dest="prefix",action="store",help='prefix of output file, default:MetaflowX',required=False,default='MetaflowX')

    return vars(p.parse_args())
    
def read_files(profile_file, taxonomy_file):
    """
    Read the input file.
    """
    profile_df = pd.read_csv(profile_file, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file, sep='\t', usecols=['BinID', 'GTDB_taxonomy'])
    return profile_df, taxonomy_df

def replace_ids_and_split(mapped_df, taxonomy_df):
    """
    Replace the profile ID with the species annotation results, and divide the species annotation information.
    """
    mapping_dict = taxonomy_df.set_index('BinID')['GTDB_taxonomy'].to_dict()

    if 'unmapped' in mapped_df['BinID'].values and show_unmap == "Y":
        mapping_dict['unmapped'] = 'unmapped;unmapped;unmapped;unmapped;unmapped;unmapped;unmapped'

    mapped_df.loc[:, 'BinID'] = mapped_df['BinID'].map(mapping_dict)
    mapped_df = mapped_df.rename(columns={'BinID': 'taxonomy'})
    
    split_taxonomy = mapped_df['taxonomy'].str.split(';', expand=True)
    mapped_df = mapped_df.drop('taxonomy', axis=1)
    taxo_mapped_df = pd.concat([split_taxonomy, mapped_df], axis=1)
    
    return taxo_mapped_df

def group_by_taxonomy_levels(profile_df):
    """
    Aggregate at each category level
    """
    taxonomy_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_list = list(profile_df.columns[len(taxonomy_levels):])
    profile_df.columns = taxonomy_levels + sample_list
    
    aggregated_dfs = {}
    for level in taxonomy_levels:
        extract_column = [level] + sample_list
        sub_pf = profile_df.loc[:, extract_column]
        aggregated_df = sub_pf.groupby(level).sum()
        aggregated_dfs[level] = aggregated_df
        
    return aggregated_dfs

def main(profile_file, taxonomy_file):
    profile_df, taxonomy_df = read_files(profile_file, taxonomy_file)
    mapped_df = replace_ids_and_split(profile_df, taxonomy_df)
    aggregated_dfs = group_by_taxonomy_levels(mapped_df)
    
    for level, df in aggregated_dfs.items():
        if show_unmap == "Y" :
            df.to_csv(f'{out_prefix}_{level}_abundance_unmap.csv')
            print(f'{out_prefix}_{level}_abundance_unmap.csv saved.')
        else:
            df.to_csv(f'{out_prefix}_{level}_abundance.csv')
            print(f'{out_prefix}_{level}_abundance.csv saved.')


if __name__ == "__main__":
    global pars,qsDir,gtdbDir,show_unmap,out_prefix
    pars = ParsReceiver()
    show_unmap = pars['show_unmap']
    out_prefix = pars['prefix']
    main(pars['profile_file'], pars['taxonomy_file'])
