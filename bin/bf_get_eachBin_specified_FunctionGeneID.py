#!/usr/bin/env python

import sys,os
import argparse as ap
from Bio import SeqIO

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    get each bin specified database functiona gene list \n
    Author:lianglifeng

    V1	2024-09-24
------------------'''
)
    ars = p.add_argument
    ars('-e',dest="function",action="store",help='Results of functional annotation using a specified database',required=True)
    ars('-d',dest="databaseName",action="store",help='pecified database name, eg VFDB, CARD....',required=True)
    ars('-i',dest="geneIDinfo",action="store",help='gene id information file',required=True)
    ars('-c',dest="geneCdhitCluster",action="store",help='gene set cd-hit cluster result file',required=True)
    ars('-b',dest="binfaFloder",action="store",help='bin fa floder',required=True)
    ars('-q',dest="outputFilePrefix",action="store",help='output file prefix default: All ',required=False,default='All')
    ars('-o',dest="outPath",action="store",help='outpath, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

# dir[allgene]= represent gene id
def getGeneCluster(filename: str) -> dict:
    f = open(filename)
    res = {}
    line = f.readline()
    while True:
        if line and line[0] == '>':
            key = None
            items = []
            for line in f:
                if not line or line[0] == '>':
                    # End of file or encounter the next cluster.
                    break
                elif '>' in line:
                    # Normal cluster content.
                    value = line.split(' ')[1][1:-3]
                    if line[-2] == '*':
                        key = value
                    items.append(value)  # Fix: Add key reactions themselves to the value list.
            if key:
                for i in items:
                    res[i] = key
            else:
                # Operations to perform when key reactions do not exist.
                pass
        else:
            # End of file.
            break
    return res


#dir[oldID] = newID
def getNewID(idF):
    with open(idF,'r') as idFile:
        for ii in idFile:
            iline = ii.strip().split('\t')
            newIDDir[iline[0]] =  iline[1]


#dir[unique]=[ko]
def pasteEggnog(eggF):
    geneFunDir={}
    with open(eggF,'r') as eggFile:
        for a in eggFile:
            if a[0] == "#":
                pass
            else:
                al = a.strip().split('\t')
                for c in  needFunction:
                    if al[c] != "-":
                        geneFunDir.setdefault(al[0],{})[c] = al[c]
                    else:
                        geneFunDir.setdefault(al[0],{})[c] = ''

    return(geneFunDir)




#dir[contig] = gene - function
def uniqueGeneFunction(eggF, idinfoF):
    uniqueGeneFunDir,contigGeneDir={},{}
    egg_geneFunctionalDir = pasteEggnog(eggF)
    with open(idinfoF,'r') as idinfoFile:
        n = 0
        for a in idinfoFile:
            if n > 0:
                al=a.strip().split('\t')
                representGene = geneClusterDir[al[-1].strip().split(' ')[0]]
                representGeneNewID = newIDDir[representGene]
                # contig - gene - function
                contigGeneDir.setdefault(al[2],{})[al[1]] = egg_geneFunctionalDir[representGeneNewID]
            n += 1


def bin2fun(binFloder,eggF,outF):
    egg_geneFunctionalDir = pasteEggnog(eggF)
    outDatabaseList=list(needFunction.keys())
    id2f={}
    for i in outDatabaseList:
        id2f[i] = open(outPath+'/'+prefix+'_'+str(needFunction[i])+'_annotation.xls', 'w')

    binFile = [os.path.join(binFloder, f) for f in os.listdir(binFloder) if f.endswith(".fa") or f.endswith(".fasta")]    #V2

    with open(outF,'w') as outFile, open(outPath+'/'+prefix+'_gene.xls','w') as geneF:
        outFile.write('binID\tnewGeneID\t'+'\t'.join([needFunction[o] for o in  outDatabaseList ])+'\n')
        for binfa in binFile:
            binfa_name = os.path.basename(binfa).rstrip(".fa")
            onebinFunDir = {}
            allbinGenelist = []

            for seq_record in SeqIO.parse(binfa, "fasta"):
                # contigID  = seq_record.id
                contigID  = "|".join(seq_record.id.split("|")[:-1]) #V2
                #print(contigID)
                if contigID in contigGeneDir: #There are contigs for which genes have not been predicted.
                    allbinGenelist = allbinGenelist + contigGeneDir[contigID]
                else:
                    print("Note : This contig does not get gene "+str(contigID))

            uniqueGene = set(filter(None, set(allbinGenelist)))
            
            print(str(binfa_name)+'\t gene number : '+str(len(uniqueGene)))
            geneF.write('%s\t%s\n'%(binfa_name,','.join(set(filter(None, set(uniqueGene))))))

            for k in uniqueGene:
                for i in needFunction:
                    if k in egg_geneFunctionalDir:
                        onebinFunDir.setdefault(i,[]).append(egg_geneFunctionalDir[k][i])
                    #else:
                    #	print("Note : Do not get this gene annatation info "+str(k))
            
            outlist=[]
            for l in outDatabaseList:
                allitemList=[]
                if l in onebinFunDir:
                    for g in onebinFunDir[l]: #bin-data:kegg
                        allitemList = allitemList + g.split(',')
                        #'\t'.join(set(filter(None,set(g.split(',').split(',')))))
                    outlist.append(','.join(set(filter(None, set(allitemList)))))
                    #print(','.join(set(filter(None, set(allitemList)))))
                    id2f[l].write('%s\t%s\n'%(binfa_name,','.join(set(filter(None, set(allitemList))))))

                else:
                    outlist.append('NA')
            outFile.write('%s\t%s\t%s\n'%(binfa_name,','.join(uniqueGene),'\t'.join(outlist))) #For each bin, only consider whether this gene exists, temporarily disregarding gene copy number.






def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outPath,prefix,shellPath,needFunction,newIDDir,geneClusterDir,contigGeneDir
    pars = ParsReceiver()
    prefix = pars['outputFilePrefix']
    makedir(pars['outPath'])
    outPath = pars['outPath']
    shellPath = os.path.split(os.path.realpath(__file__))[0]
    newIDDir = {}
    
    needFunction={1: pars['databaseName'] }


    #1 get new old gene id pair
    getNewID(pars['geneIDinfo'])

    #2 get gene cluster reprsent info
    geneClusterDir = getGeneCluster(pars['geneCdhitCluster'])

    allgeneNewDir,contigGeneDir={},{}
    for g in geneClusterDir:
        re = geneClusterDir[g]
        re_newID = newIDDir.get(re, None)  # 安全访问
        if re_newID is None:
            continue  # 跳过或记录日志
        
        allgeneNewDir[g] = re_newID
        gg =  g.strip().split("|")
        contigID ="%s|%s" % (gg[0],gg[-1].split('_')[0])
        contigGeneDir.setdefault(contigID,[]).append(re_newID)

    
    #3 get contig gene info
    #uniqueGeneFunction(pars['eggnog'],pars['geneIDinfo'])

    #4 out put each bin function
    bin2fun(pars['binfaFloder'],pars['function'],outPath+'/'+prefix+'.xls')


if __name__ == '__main__':
    main()