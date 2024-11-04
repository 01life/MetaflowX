#!/usr/bin/env python

import sys,re,os
import argparse as ap

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    Funtion:
        stat each bin different database function catagory to plot
    Version:
        V1	2023-07-07
    Author:
        Liang Lifeng
------------------a'''
)
    ars = p.add_argument
    ars('-k',dest="ko",action="store",help='eachbin*KEGG_ko.txt',required=True)
    ars('-K',dest="KEGG",action="store",help='ko_category.txt',required=True)

    ars('-g',dest="go",action="store",help='eachbin*GOs.txt',required=True)
    ars('-G',dest="GO",action="store",help='GO.ontology.txt',required=True)

    ars('-c',dest="cog",action="store",help='eachbin*cog_catF.txt',required=True)
    ars('-C',dest="COG",action="store",help='cog.level.txt',required=True)

    ars('-a',dest="cazy",action="store",help='eachbin*CAZy.txt',required=True)
    ars('-A',dest="CAZY",action="store",help='cazy_category.txt',required=True)
    return vars(p.parse_args())

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)



def patseDataBase(profleF, CategoryF, outF, databseName):
    #outF: binID \t database \t category \t count \n
    oneCategorieDir, oneDir={},{}
    with open(profleF,'r') as oneProfileF, open(CategoryF,'r') as oneDBcategoryF:
        for o in oneDBcategoryF:
            ol = o.strip().split('\t')
            oneCategorieDir[ol[0]] = ol[-1]
        for i in oneProfileF:
            iline = i.strip().split('\t')
            oneDir = {}
            #eachBin[each category] = count
            if len(iline) >1:
                for b in iline[1].split(','):
                    if databseName == 'COG':
                        for d in b:
                            if d not in oneCategorieDir:
                                oneDir.setdefault('Others',[]).append(d)
                            else:
                                oneDir.setdefault(oneCategorieDir[d],[]).append(d)

                    else:
                        if b not in oneCategorieDir:
                            oneDir.setdefault('Others',[]).append(b)
                        else:
                            oneDir.setdefault(oneCategorieDir[b],[]).append(b)

                for e in oneDir:
                    outF.write("%s\t%s\t%s\t%s\n"%(iline[0],databseName,str(e),str(len(oneDir[e]))))
            else:
                print("Waring~~~~~")
                print(str(iline[0])+" bin didn't get "+str(databseName)+' result.')
 
def pasteCAZY(cazyF,CategoryF,outF):

    oneDir, cazyCategorieDir ={},{}
    with open(cazyF,'r') as oneProfileF, open(CategoryF,'r') as CategoryFile:
        for c in CategoryFile:
            cl = c.strip().split("\t")
            cazyCategorieDir[cl[0]] = cl[1]

        for i in oneProfileF:
            iline = i.strip().split('\t')
            oneDir = {}
            #eachBin[each category] = count
            if len(iline) >1:
                for b in iline[1].split(','):
                    cazyClass = b[:2]
                    if cazyClass not in cazyCategorieDir:
                        oneDir.setdefault('Others',[]).append(cazyClass)
                    else:
                        oneDir.setdefault(cazyCategorieDir[cazyClass],[]).append(cazyClass)

            for e in oneDir:
                outF.write("%s\t%s\t%s\t%s\n"%(iline[0],'CAzY',str(e),str(len(oneDir[e]))))


def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outpath,outFile
    pars = ParsReceiver()

    outFile = open('./eachBin_Function_stat.xls','w')
    outFile.write('binID\tdatabase\tclass\tcount\n')
    patseDataBase(pars['ko'],pars['KEGG'] , outFile, 'KEGG')
    patseDataBase(pars['go'], pars['GO'], outFile, 'GO')
    patseDataBase(pars['cog'], pars['COG'], outFile, 'COG')
    pasteCAZY(pars['cazy'],pars['CAZY'],outFile)
    outFile.close()


if __name__ == '__main__':
    main()

