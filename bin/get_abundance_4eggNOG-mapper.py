#!/usr/bin/env python

import sys,re,os
import argparse as ap


def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    get special function database profile from eggNOG-mapper result \n
    Author:lianglifeng

    V1	2023-03-20
    V2  2024-09-04 check nagative result
------------------'''
)
    ars = p.add_argument
    ars('-e',dest="eggnog",action="store",help='*emapper.annotations from eggNOG-mapper',required=True)
    ars('-f',dest="profile",action="store",help='gene abundance profile',required=True)
    ars('-g',dest="genePrefix",action="store",help='prefix of geneID, letter only default:NC ',required=False,default='NC')
    ars('-q',dest="outputFilePrefix",action="store",help='output file prefix default: All ',required=False,default='All')
    ars('-o',dest="outPath",action="store",help='outpath, default:./',required=False,default=os.getcwd())
    return vars(p.parse_args())

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def check_file(filename):
    # Check if the number of lines in the file is greater than 2.
    with open(filename, 'r') as file:
        lines = file.readlines()
        if len(lines) > 1:
            return True
        else:
            return False


def dumpWrite(geneid, className):
    dd = ''
    for d in className:
        dd = dd + str(geneid) + '\t' + str(d) + '\n'
    return(dd)


def reNameProfile(geneProfie,outPath,prefix,genePrefix):
    newF = open(outPath+'/'+prefix+'_geneset_gene_abundance.xls','w')
    with open(geneProfie,'r') as oldF:
        n = 0
        for g in oldF:
            if n ==0:
                newF.write(str(g))
            else:
                gl=g.strip().split('\t')
                newID  = str(genePrefix)+"_{0:010d}".format(int(gl[0]))
                newF.write(str(newID)+'\t'+'\t'.join(gl[1:])+'\n')
            n +=1
    newF.close()
    
            

def pasteEggnog(eggF,profileF,outPath,prefix):
    needFunction={6:"cog_catF",9:"GOs",10:"EC",11:"KEGG_ko",12:"KEGG_Pathway",13:"KEGG_Module", 14:"KEGG_Reaction",18:"CAZy",20:"PFAMs"}
    id2f = {}
    for i in needFunction:
        id2f[i] = open(outPath+'/'+prefix+'_geneset_function_'+str(needFunction[i])+'_annotation.xls', 'w')
    
    with open(eggF,'r') as eggFile:
        for a in eggFile:
            if a[0] == '#':
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
                            cc = dumpWrite(al[0], newC)
                        elif c == 12:
                            orgc = al[c].strip().split(',')
                            newC = []
                            for d in orgc:
                                if re.search('^map',d):
                                    newC.append(d)
                            cc = dumpWrite(al[0], newC)
                        else:
                            cc = dumpWrite(al[0], al[c].strip().split(','))
                        id2f[c].write(cc)
    
    for d in  needFunction:
        id2f[d].close()
        filename = os.path.abspath(outPath+'/'+prefix+'_geneset_function_'+str(needFunction[d])+'_annotation.xls')
        if check_file(filename):
            commond = "perl %s/functional_profile.pl -k %s -f %s -o %s -pr %s_geneset_function_%s_abundance.xls " %(shellPath,filename,os.path.abspath(profileF),outPath,prefix,needFunction[d])
            os.system(commond)
            print("Runing "+ needFunction[d]+'  ~~~')

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outPath,prefix,shellPath,genePrefix
    pars = ParsReceiver()
    prefix = pars['outputFilePrefix']
    makedir(pars['outPath'])
    outPath = pars['outPath']
    
    shellPath = os.path.split(os.path.realpath(__file__))[0]
    genePrefix = pars['genePrefix']
    reNameProfile(pars['profile'],outPath,prefix,genePrefix)
    pasteEggnog(pars['eggnog'],outPath+'/'+prefix+'_geneset_gene_abundance.xls',outPath,prefix)

if __name__ == '__main__':
    main()