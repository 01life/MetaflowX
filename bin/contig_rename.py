#!/usr/bin/env python

import sys,re,os
import argparse as ap
from Bio import SeqIO

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    rename contig and return new-old ID relationship \n
    Version:
        V1	2023-09-04
    Author:
        Lifeng Liang
------------------'''
)
    ars = p.add_argument
    ars('-t',dest="oldfa",action="store",help='old contig fa',required=True)
    ars('-s',dest="sampleID",action="store",help='sample ID ',required=True)
    ars('-m',dest="mapping",action="store",help='id relationship mapping outfile ',required=True)
    ars('-o',dest="outFile",action="store",help='new fa file ',required=True)
    return vars(p.parse_args())


def renameSeq(contigFa,sampleID,outFa,mapF):
    n = 1
    with open(outFa,'w')  as outFafile, open(mapF,'w') as mapFile:
        for seq_record in SeqIO.parse(contigFa, "fasta"):
            newID = sampleID + "|"+str(n)
            outFafile.write(f">{newID}\n{seq_record.seq}\n")
            mapFile.write(f"{seq_record.id}\t{newID}\n")
            n += 1


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outpath,contigLenDir
    pars = ParsReceiver()
    contigLenDir = {}
    renameSeq(pars['oldfa'],pars['sampleID'],pars['outFile'],pars['mapping'])


if __name__ == '__main__':
    main()
