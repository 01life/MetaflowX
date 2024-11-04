#!/usr/bin/env python

import sys,re,os
import argparse as ap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    filter gene length  & rename geneset cds ID & give the cds info\n
    Author:lianglifeng

    V1	2023-03-17
------------------'''
)
    ars = p.add_argument
    ars('-n',dest="nucleotide",action="store",help='nucleotide sequences fa file ',required=True)
    ars('-p',dest="protein",action="store",help='protein  sequences fa file  ',required=True)
    ars('-g',dest="genePrefix",action="store",help='prefix of geneID, letter only default:NC ',required=False,default='NC')
    ars('-q',dest="outputFilePrefix",action="store",help='output file prefix default: All ',required=False,default='All')
    ars('-o',dest="outPath",action="store",help='outpath, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def renameID(cdsF,pepF):
    with open(outPath+'/'+prefix+'_geneset_gene_info.xls','w') as allgeneInfoF, open(outPath+'/'+prefix+'_geneset_gene_length.xls','w') as filtergeneInfoF, open (outPath+'/'+prefix+'_geneset_gene.fa','w') as outcdsF, open( outPath+'/'+prefix+'_geneset_protein.fa','w') as outpepF:
        allgeneInfoF.write('Number\tnewID\tcontigID\toriginalDescription\n')
        filtergeneInfoF.write('ID\tname\tlength\n')
        
        need_record,n,idDir,idlist = [],1,{},set()
        #handle cds
        for seq_record in SeqIO.parse(cdsF, "fasta"):
            contigID = seq_record.description.split("#")[0].strip()
            rawID =  seq_record.description
            newID = str(genePrefix)+"_{0:010d}".format(n)
            seqLen = len(seq_record.seq)
            idDir[rawID] = newID
            idlist.add(rawID)
            allgeneInfoF.write('\t'.join([str(n),newID,contigID,rawID])+'\n')
            filtergeneInfoF.write('\t'.join([str(n),newID,str(seqLen)])+'\n')
            outcdsF.write('>' +str(newID)+'\n'+str(seq_record.seq)+'\n')
            n += 1

        #handle pep
        needpep=[]
        for pep in SeqIO.parse(pepF, "fasta"):
            pepID = pep.description
            if pepID in idlist:
                outpepF.write('>' +str(idDir[pepID])+'\n'+str(pep.seq)+'\n')


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars, genePrefix,outPath,prefix
    pars = ParsReceiver()
    genePrefix = pars['genePrefix']
    prefix = pars['outputFilePrefix']
    makedir(pars['outPath'])
    outPath = pars['outPath']
    renameID(pars['nucleotide'],pars['protein'])

if __name__ == '__main__':
    main()
