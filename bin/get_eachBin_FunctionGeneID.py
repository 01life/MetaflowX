#!/usr/bin/env python

import sys,re,os
import argparse as ap
from Bio import SeqIO

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    get each bin functiona gene list \n
    Author:lianglifeng xiaoyang

    V1	2023-03-21
    V2  2024-09-30 using os
------------------'''
)
    ars = p.add_argument
    ars('-e',dest="eggnog",action="store",help='*emapper.annotations from eggNOG-mapper',required=True)
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
                    break
                elif '>' in line:
                    value = line.split(' ')[1][1:-3]
                    if line[-2] == '*':
                        key = value
                    items.append(value)  
            if key:
                for i in items:
                    res[i] = key
            else:
                pass
        else:
            break
    return res


#dir[oldID] = newID
def getNewID(idF):
    with open(idF,'r') as idFile:
        for ii in idFile:
            iline = ii.strip().split('\t')
            newIDDir[iline[-1]] =  iline[0]


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
                        if c == 11:
                            orgc = al[c].strip().split(',')
                            newC = []
                            for d in orgc:
                                newC.append(d.strip().split(":")[-1])
                            geneFunDir.setdefault(al[0],{})[c] = ','.join(newC)
                        elif c == 12:
                            orgc = al[c].strip().split(',')
                            newC = []
                            for d in orgc:
                                if re.search('^map',d):
                                    newC.append(d)
                            geneFunDir.setdefault(al[0],{})[c] = ','.join(newC)
                        else:
                            newC = al[c].strip().split(',')
                            geneFunDir.setdefault(al[0],{})[c] = ','.join(newC)
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
                contigID  = "|".join(seq_record.id.split("|")[:-1]) #V2

                if contigID in contigGeneDir:
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

            
            outlist=[]
            for l in outDatabaseList:
                allitemList=[]
                if l in onebinFunDir:
                    for g in onebinFunDir[l]: #bin-data:kegg
                        allitemList = allitemList + g.split(',')
                    outlist.append(','.join(set(filter(None, set(allitemList)))))
                    id2f[l].write('%s\t%s\n'%(binfa_name,','.join(set(filter(None, set(allitemList))))))

                else:
                    outlist.append('NA')
            outFile.write('%s\t%s\t%s\n'%(binfa_name,','.join(uniqueGene),'\t'.join(outlist)))

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
    
    needFunction={6:"cog_catF",9:"GOs",10:"EC",11:"KEGG_ko",12:"KEGG_Pathway",13:"KEGG_Module", 14:"KEGG_Reaction",18:"CAZy",20:"PFAMs"}

    commond1 = 'cut -f 2,3 ' + pars['geneIDinfo'] + '> '+ pars['geneIDinfo']+'.tmp'
    os.system(commond1)

    #1 get new old gene id pair
    getNewID(pars['geneIDinfo']+'.tmp')

    #2 get gene cluster reprsent info
    geneClusterDir = getGeneCluster(pars['geneCdhitCluster'])

    allgeneNewDir,contigGeneDir={},{}
    for g in geneClusterDir:
        re=geneClusterDir[g]
        re_newID = newIDDir[re]
        allgeneNewDir[g] = re_newID
        gg =  g.strip().split("|")
        contigID ="%s|%s" % (gg[0],gg[-1].split('_')[0])
        contigGeneDir.setdefault(contigID,[]).append(re_newID)

    
    #3 get contig gene info
    #uniqueGeneFunction(pars['eggnog'],pars['geneIDinfo'])

    #4 out put each bin function
    bin2fun(pars['binfaFloder'],pars['eggnog'],outPath+'/'+prefix+'.xls')


if __name__ == '__main__':
    main()