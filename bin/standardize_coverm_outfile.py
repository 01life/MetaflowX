#!/usr/bin/env python

import sys
import argparse as ap
import shutil
from pathlib import Path

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    Funtion:
        standardize coverm resule file column headers \n
    Version:
        V1	2024-04-28
    Author:
        Lifeng Liang
------------------'''
)
    ars = p.add_argument
    ars('-i',dest="input",action="store",help='input file',required=True)
    ars('-o',dest="outfile",action="store",help='outfile',required=True)
    ars('-s',dest="sortedBam",action="store",help='using sorted bam to calculate the abundace or Not.(Y/N)[Default=Y]',required=False,default='Y')
    ars('-m',dest="method",action="store",help='CoverM calculation method, needs to be similar to the CoverM option.[mean,relative_abundance,trimmed_mean,covered_bases,variance,length,count,rpkm,tpm,reads_per_base]',required=True)
    return vars(p.parse_args())


def renameF(inputFile, outputFile, method, sortedbam):
    with open(inputFile,'r') as inputF, open(outputFile,'w') as outputF:
        
        if sortedbam == "Y":
            resubName = '.sorted '+str(renameDir[method])
        else:
            resubName = str(renameDir[method])

        header = inputF.readline().strip().split("\t")
        newheader = []
        for h in header[1:]:
            newID = h.strip().split(resubName)[0].strip()
            newheader.append(newID)
        outputF.write("BinID\t" + "\t".join(newheader)+'\n')

        for  a in inputF:
            outputF.write(a)


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,renameDir
    pars = ParsReceiver()
    
    renameDir = {
        "mean":"Mean",
        "trimmed_mean":"Trimmed Mean",
        "covered_bases":"Covered Bases",
        "rpkm":"RPKM",
        "tpm":"TPM",
        "variance":"Variance",
        "length":"Length",
        "count":"Read Count",
        "reads_per_base":"Reads per base",
        "relative_abundance":"Relative Abundance (%)",
        }
    

    if pars['method'] not in renameDir:
        print(f" -m method name do not in CoverM method, Please check! [mean,relative_abundance,trimmed_mean,covered_bases,variance,length,count,rpkm,tpm,reads_per_base] ")
        print(f"only change the file name! ")

        shutil.copy(Path(pars['input']), Path(pars['outfile']))

    else:
        renameF(pars['input'], pars['outfile'],pars['method'],pars['sortedBam'])


if __name__ == '__main__':
    main()