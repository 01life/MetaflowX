#!/usr/bin/env python

import sys,Bio
from Bio.Seq import Seq
from Bio import SeqIO

if float(Bio.__version__) > 1.8:
    from Bio.SeqUtils import gc_fraction as gc
else:
    from Bio.SeqUtils import GC as gc


with open(sys.argv[2],'w') as outF:
    outF.write('ContigID\tSample\tLength\tGC\n')
    for record in SeqIO.parse(sys.argv[1], 'fasta'):
        sampleID = record.id.split("|")[0]
        seqlen = len(record.seq)
        contig_gc = '%.4f'%gc(record.seq)
        outF.write(str(record.id)+"\t"+str(sampleID)+'\t'+str(seqlen)+'\t'+str(contig_gc)+'\n')

