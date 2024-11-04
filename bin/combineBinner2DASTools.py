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



def parseContig2bin(tsvF,outF,contigDir,pepDir):
    sampleList,binnerList =set(),set()
    sampleBinnerDir={}
    combinations = []
    with open(tsvF,'r') as tsvFile, open(outF,'w') as outFile:
        for a in tsvFile:
            sampleID, binnerName, tsvPath = a.strip().split("\t")
            sampleList.add(sampleID)
            binnerList.add(binnerName)
            sampleBinnerDir.setdefault(sampleID,{})[binnerName] = os.path.abspath(tsvPath)
        
        if len(binnerList) >= 2:
            for r in range(2, len(binnerList) + 1):
                combinations.extend(list(itertools.combinations(binnerList, r)))
        else:
            print("binner has less than 2 trees, DASTools not applicable. Program terminated.")
            sys.exit()
            

        for sample in sampleList:
            for oneCombination in combinations:
                tsvlist = [sampleBinnerDir[sample][binner] for binner in oneCombination]
                combinetxt = ",".join(tsvlist)
                binnertxt = "_".join(oneCombination)
                if sample in contigDir and sample in pepDir:
                    outFile.write(f"{sample}\t{contigDir[sample]}\t{pepDir[sample]}\t{binnertxt}\t{combinetxt}\n")
                else:
                    print(f"do not get the protein file or contigfile of sample: {sample}, please check!!")
                    sys.exit()
                


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

    parseContig2bin(pars['contig2bin'],outfile,contigDir,pepDir)

if __name__ == '__main__':
    main()

