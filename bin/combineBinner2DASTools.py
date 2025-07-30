#!/usr/bin/env python

import sys,os,itertools
import argparse as ap

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    combine multi binner result input to DASTools\n
    Version:
        V1	2023-12-20
    Author:
        Lifeng Liang
------------------'''
)
    ars = p.add_argument
    ars('-t',dest="contig2bin",action="store",help='sample binner contig2bin.tsv list',required=True)
    ars('-p',dest="pep",action="store",help='sample protein sequence fa list',required=True)
    ars('-c',dest="contig",action="store",help='sample contig sequence fa list',required=True)
    ars('-o',dest="outPath",action="store",help='outpathh, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def getContig(contigF):
    with open(contigF,'r') as contigFile:
        for c in contigFile:
            sampleID, contigPath = c.strip().split('\t')
            contigDir[sampleID] = contigPath

def getPep(pepF):
    with open(pepF,'r') as pepFile:
        for c in pepFile:
            sampleID, pepPath = c.strip().split('\t')
            pepDir[sampleID] = pepPath


def parse_contig2bin(tsv_path: str, output_file: str, contig_map: dict, pep_map: dict):
    """Generate binner combination TSV file for DASTools."""
    sample_binner_map = {}
    sample_ids = set()
    binner_names = set()

    with open(tsv_path, 'r') as f:
        for line in f:
            sample_id, binner, bin_path = line.strip().split('\t')
            sample_ids.add(sample_id)
            binner_names.add(binner)
            sample_binner_map.setdefault(sample_id, {})[binner] = os.path.abspath(bin_path)

    if len(binner_names) < 2:
        print("Less than 2 binners provided. DASTools requires multiple binners. Exiting.")
        sys.exit(1)

    combinations = []
    for r in range(2, len(binner_names) + 1):
        combinations.extend(itertools.combinations(binner_names, r))

    with open(output_file, 'w') as out_f:
        for sample_id in sample_ids:
            if sample_id not in contig_map or sample_id not in pep_map:
                print(f"[WARNING] Missing contig or protein file for sample: {sample_id}. Skipped.")
                continue

            for comb in combinations:
                if not all(binner in sample_binner_map[sample_id] for binner in comb):
                    missing = [b for b in comb if b not in sample_binner_map[sample_id]]
                    # print(f"[INFO] Skipping combination {comb} for sample {sample_id} due to missing binner(s): {', '.join(missing)}")
                    continue

                bin_paths = [sample_binner_map[sample_id][b] for b in comb]
                combined_bins = ",".join(bin_paths)
                combo_name = "_".join(comb)
                out_f.write(f"{sample_id}\t{contig_map[sample_id]}\t{pep_map[sample_id]}\t{combo_name}\t{combined_bins}\n")


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outpath,contigDir,pepDir
    pars = ParsReceiver()
    contigDir,pepDir = {},{}
    makedir(pars['outPath'])
    outpath = pars['outPath']


    getContig(pars['contig'])
    getPep(pars['pep'])


    outfile = f"{outpath}/binner_combination.txt"

    parse_contig2bin(pars['contig2bin'],outfile,contigDir,pepDir)

if __name__ == '__main__':
    main()

