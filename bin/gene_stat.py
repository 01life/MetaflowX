#!/usr/bin/env python

from Bio import SeqIO
import sys

sampleGeneDir={}
with open(sys.argv[2],'w') as outF:
    outF.write('Sample\tgeneNumber\tgeneAverageLength\n')
    for record in SeqIO.parse(sys.argv[1], 'fasta'):
        sampleID = record.id.split("|")[0]
        sampleGeneDir.setdefault(sampleID,[]).append(len(record.seq))

    for a in sampleGeneDir.keys():
        valuelist = sampleGeneDir[a]
        geneNumber = len(valuelist)
        avgGeneLength = sum(valuelist) / geneNumber
        outF.write("%s\t%s\t%s\n"%(a,str(geneNumber),str('%.4f'%(avgGeneLength))))