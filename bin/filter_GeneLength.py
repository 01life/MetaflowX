#!/usr/bin/env python

import sys
import os
import argparse as ap
from Bio import SeqIO
import math

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    filter gene length\n
    Author:lianglifeng

    V1	2023-03-17
------------------'''
)
    ars = p.add_argument
    ars('-i',dest="inputFa",action="store",help='raw fa file ',required=True)
    ars('-s',dest="splitThreshold",action="store",help='the threshold of split run cd-hit, default[8000000]',required=False,default=8000000)
    ars('-c',dest="chunkSize",action="store",help='gene count of each cd-hit split task, default[300000]',required=False,default=300000)
    ars('-l',dest="mimlen",action="store",help='minimum length of gene default:150',required=False,default='150')
    ars('-o',dest="outFa",action="store",help=' output fa file',required=True)
    return vars(p.parse_args())


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def fillterLength(inputFa, outFa, mimlen,splitThreshold,chunkSize,):
    with open(inputFa+'.all.geneLength.txt','w') as infoF, open(outFa,'w') as outFaFile:
        infoF.write('ID\tLength\n')
        need_record,minLen, n= [],float(mimlen),0
        for seq_record in SeqIO.parse(inputFa, "fasta"):
            sl = len(seq_record.seq)
            infoF.write(str(seq_record.id)+'\t'+str(sl)+'\n')
            if sl >= minLen:
                n+=1
                outFaFile.write(">"+seq_record.description+'\n'+str(seq_record.seq)+'\n')

        if n > int(splitThreshold) :

            splitNum = math.ceil(n / int(chunkSize))
            splitFile = open(str(splitNum)+'.multi.cdhit.task','w')
            splitFile.close()

        else:
            splitFile = open('single.cdhit.task.num','w')
            splitFile.close()
        
def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,mimlen,shellPath,queueType
    pars = ParsReceiver()
    fillterLength(pars['inputFa'],pars['outFa'],pars['mimlen'],pars['splitThreshold'],pars['chunkSize'])

if __name__ == '__main__':
    main()
