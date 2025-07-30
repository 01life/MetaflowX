#!/usr/bin/env python

import sys,os,itertools
import argparse as ap

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    create multi binner result input to DASTools\n
    Version:
        V1	2025-07-22
    Author:
        Lifeng Liang
------------------'''
)
    ars = p.add_argument
    ars('-t',dest="contig2bin",action="store",help='sample binner contig2bin.tsv list',required=True)
    ars('-o',dest="outPath",action="store",help='outpathh, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def parse_contig2bin(tsv_path: str, outdir: str):
    sample_binner_map = {}

    with open(tsv_path, 'r') as f:
        for line in f:
            sample_id, binner, bin_path = line.strip().split(',')
            sample_binner_map.setdefault(sample_id, {})[binner] = os.path.abspath(bin_path)

    for sample_id, binner_bins in sample_binner_map.items():
        sample_outdir = os.path.join(outdir, sample_id)
        makedir(sample_outdir)

        # link  contig2bin 
        for binner, bin_file in binner_bins.items():
            if not os.path.exists(bin_file):
                print(f"[WARNING] Bin file not found: {bin_file}. Skipped.")
                continue

            link_name = os.path.join(sample_outdir, f"{sample_id}_{binner}.tsv")
            try:
                if not os.path.exists(link_name):
                    os.symlink(bin_file, link_name)
            except Exception as e:
                print(f"[ERROR] Linking bin file failed: {bin_file} -> {link_name}: {e}")



def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outpath
    pars = ParsReceiver()
    makedir(pars['outPath'])
    outpath = pars['outPath']
    parse_contig2bin(pars['contig2bin'],outpath)

if __name__ == '__main__':
    main()

