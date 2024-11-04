#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
import sys


with open(sys.argv[1],'r') as inputF, open(sys.argv[2],'w') as outF:
    outF.write("BinID\tcontigNumber\tGenomeSize\n")
    for a in inputF:
        al = a.strip().split('\t')
        contigNum, genomeSize = 0,0
        for record in SeqIO.parse(al[1], 'fasta'):
            contigNum += 1
            genomeSize = genomeSize + len(record.seq)
        outF.write(str(al[0])+'\t'+str(contigNum)+'\t'+str(genomeSize)+'\n')
