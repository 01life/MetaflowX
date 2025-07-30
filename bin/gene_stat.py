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
    V2  2025-02-26
------------------'''
)
    ars = p.add_argument
    ars('-n',dest="nucleotide",action="store",help='nucleotide sequences fa file ',required=True)
    ars('-p',dest="protein",action="store",help='protein  sequences fa file  ',required=True)
    ars('-l',dest="mimlen",action="store",help='minimum length of gene default:150',required=False,default='150')
    ars('-g',dest="genePrefix",action="store",help='prefix of geneID, letter only default:NC ',required=False,default='NC')
    ars('-q',dest="outputFilePrefix",action="store",help='output file prefix default: sampleID ',required=False,default='All')
    ars('-o',dest="outPath",action="store",help='outpath, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def renameID(cdsF,pepF):
    with open(outPath+'/'+prefix+"_L"+minilen+'_gene_info.xls','w') as allgeneInfoF, open(outPath+'/'+prefix+"_L"+minilen+'_gene_length.xls','w') as filtergeneInfoF, open (outPath+'/'+prefix+"_L"+minilen+'_cds.fa','w') as outcdsF, open( outPath+'/'+prefix+"_L"+minilen+'_protein.fa','w') as outpepF, open(outPath+'/'+prefix+"_L"+minilen+'_gene_stat.xls','w') as allgeneStatF: 
        allgeneInfoF.write('Number\tgeneID\tcontigID\toriginalDescription\n')
        filtergeneInfoF.write('ID\tname\tlength\n')
        
        need_record,n,idDir,idlist = [],1,{},set()
        newID_cds = {}
        newID_cds_len = {}
        #handle cds
        for seq_record in SeqIO.parse(cdsF, "fasta"):
            seqLen = len(seq_record.seq)
            if seqLen > int(minilen):
                contigID = seq_record.description.split("#")[0].strip()
                rawID =  seq_record.description
                newID = seq_record.id
                # newID = str(genePrefix)+"_{0:010d}".format(n)
                
                idDir[rawID] = newID
                idlist.add(rawID)
                newID_cds_len[newID] = seqLen
                newID_cds[newID] = seq_record.seq

                # allgeneInfoF.write('\t'.join([str(n),newID,contigID,rawID])+'\n')
                # filtergeneInfoF.write('\t'.join([str(n),newID,str(seqLen)])+'\n')
                allgeneInfoF.write('\t'.join([str(n),newID,contigID,rawID])+'\n')
                filtergeneInfoF.write('\t'.join([str(n),newID,str(seqLen)])+'\n')

                n += 1

        #handle pep
        newID_pep_seq={}
        for pep in SeqIO.parse(pepF, "fasta"):
            pepID = pep.description
            if pepID in idlist:
                newID_pep_seq[idDir[pepID]] = pep.seq
                

        #sort by length
        newID_cds_len_sort = sorted(newID_cds_len.items(), key=lambda x:x[1], reverse=True)
        for newID,seqLen in newID_cds_len_sort:
            outcdsF.write('>' +str(newID)+'\n'+str(newID_cds[newID])+'\n')
            outpepF.write('>' +str(newID)+'\n'+str(newID_pep_seq[newID])+'\n')

        # stat the whole gene set
        allgeneStatF.write('Sample\tgeneNumber\tgeneAverageLength\n')

        valuelist = newID_cds_len.values()
        geneNumber = len(newID_cds_len.keys())
        avgGeneLength = sum(valuelist) / geneNumber
        allgeneStatF.write("%s\t%s\t%s\n"%(prefix,str(geneNumber),str('%.4f'%(avgGeneLength))))

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars, genePrefix,outPath,prefix,minilen
    pars = ParsReceiver()
    genePrefix = pars['genePrefix']
    prefix = pars['outputFilePrefix']
    makedir(pars['outPath'])
    outPath = pars['outPath']
    minilen = pars['mimlen']
    renameID(pars['nucleotide'],pars['protein'])

if __name__ == '__main__':
    main()
