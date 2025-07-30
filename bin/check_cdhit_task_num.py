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
    check cdhit task number\n
    Author:lianglifeng

    V1	2023-03-17
    V2  2025-02-26
------------------'''
)
    ars = p.add_argument
    ars('-i',dest="stat",action="store",help='all sample gene stat info',required=True)
    ars('-s',dest="splitThreshold",action="store",help='the threshold of split run cd-hit, default[8000000]',required=False,default=8000000)
    ars('-c',dest="chunkSize",action="store",help='gene count of each cd-hit split task, default[300000]',required=False,default=300000)
    return vars(p.parse_args())


def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def split_task(statF,splitThreshold,chunkSize,):

    genenumber = 0
    with open(statF,'r') as statFile:
        header = statFile.readline()
        for line in statFile:
            allinfo = line.strip().split('\t')
            genenumber  += int(allinfo[1])


        if genenumber > int(splitThreshold) :

            splitNum = math.ceil(genenumber / int(chunkSize))
            splitFile = open(str(splitNum)+'.multi.cdhit.task','w')
            splitFile.close()

        else:
            splitFile = open('single.cdhit.task','w')
            splitFile.close()
        
def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,shellPath,queueType
    pars = ParsReceiver()
    split_task(pars['stat'],pars['splitThreshold'],pars['chunkSize'])

if __name__ == '__main__':
    main()
