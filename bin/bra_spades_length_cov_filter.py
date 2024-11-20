#!/usr/bin/env python

import sys
import argparse as ap
from Bio import SeqIO
import numpy as np
import pandas as pd


def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    base on the spades orginal contig ID filter contig\n
    Version:
        V1	2024-07-02
    Author:
        Lifeng Liang
------------------'''
)
    ars = p.add_argument
    ars('-f',dest="oldfa",action="store",help='old contig fa',required=True)
    ars('-l',dest="len",action="store",help='minimum length of the contig ',required=False,default=2000,type=int)
    ars('-p',dest="prefix",action="store",help='output prefix ',required=True)
    return vars(p.parse_args())

def get_contig_info(contigFa,infoFile,minilength):
    with open(infoFile,'w')  as outF:
        outF.write('seq.id\tcov\tlength\n')
        for seq_record in SeqIO.parse(contigFa, "fasta"):
            #NODE_1_length_362809_cov_5.671360
            id_info = seq_record.id.strip().split("_")
            contig_cov = id_info[-1]
            contig_len = int(id_info[-3])
            if contig_len >= minilength:
                outF.write(f"{seq_record.id}\t{contig_cov}\t{contig_len}\n")


def log_weighted_zscore_outliers(data, weights, threshold=3, epsilon=1e-10):
    log_data = np.log(data + epsilon) # Avoid log(0)
    weighted_mean = np.average(log_data, weights=weights)
    weighted_var = np.average((log_data - weighted_mean)**2, weights=weights)
    weighted_std = np.sqrt(weighted_var)
    
    # Avoid division by zero
    if weighted_std == 0:
        return np.zeros(len(data), dtype=bool)
    
    weighted_z_scores = (log_data - weighted_mean) / weighted_std
    return np.abs(weighted_z_scores) > threshold

def modified_zscore_outliers(data, threshold=3.5, epsilon=1e-10):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    
    # Avoid division by zero
    if mad == 0:
        return np.zeros(len(data), dtype=bool)
    
    modified_z_scores = 0.6745 * (data - median) / mad
    return np.abs(modified_z_scores) > threshold

def iqr_outliers(data, factor=1.5):
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - factor * IQR
    upper_bound = Q3 + factor * IQR

    return (data < lower_bound) | (data > upper_bound)


def remove_coverage_outlier_by_weighted(input_file, z_score_threshold=3):
    # Read data.
    contig_coverage_df = pd.read_csv(input_file, sep="\t", index_col=0)
    
    #IQR method.
    iqr_out = iqr_outliers(contig_coverage_df['cov'])
    #Modified Z-score method.
    mod_zscore_out = modified_zscore_outliers(contig_coverage_df['cov'])

    #Weighted Z-score method after log transformation.
    log_weighted_zscore_out = log_weighted_zscore_outliers(contig_coverage_df['cov'], weights=contig_coverage_df['length'], threshold=z_score_threshold)

    # Merge the results of all methods.

    print(f"Weighted Z-score method found {sum(log_weighted_zscore_out)} outliers")
    print(f"IQR method found {sum(iqr_out)} outliers")
    print(f"Modified Z-score method found {sum(mod_zscore_out)} outliers")

    # all_outliers =  iqr_out | log_weighted_zscore_out | mod_zscore_out
    all_outliers = iqr_out & log_weighted_zscore_out & mod_zscore_out

    print(f"mergar all method {sum(all_outliers)} outliers")

    # Remove all outliers from contig_coverage_df.
    contig_coverage_df_filtered = contig_coverage_df[~all_outliers]

    return contig_coverage_df_filtered.index.tolist()

def get_out_fa(oldFa,keep_idlist,outFa):
    with open(outFa,'w')  as outFaFile:
        for seq_record in SeqIO.parse(oldFa, "fasta"):
            if seq_record.id in keep_idlist:
                outFaFile.write(f'>{seq_record.id}\n{seq_record.seq}\n')


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars
    pars = ParsReceiver()
    get_contig_info(pars['oldfa'],pars['prefix']+'_contig_info.txt',pars['len'])
    keep_seq_id = remove_coverage_outlier_by_weighted(pars['prefix']+'_contig_info.txt')
    get_out_fa(pars['oldfa'],keep_seq_id,pars['prefix']+'_remove_outlier.fa')


if __name__ == '__main__':
    main()
