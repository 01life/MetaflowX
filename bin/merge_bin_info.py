#!/usr/bin/env python
import sys
import argparse as ap

def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    
    Funtion:
        merge bin information to show in report\n
    Version:
        V1	2024-03-08
    Author:
        Lifeng Liang
-------------------'''
)
    ars = p.add_argument
    ars('-m',dest="map",action="store",help='MetaFlowX_HQ_unique_bins_rename_map.xls',required=True)
    ars('-g',dest="gtdb",action="store",help='gtdbtk.taxonomy2ncbi.summary.tsv',required=True)
    ars('-s',dest="summary",action="store",help='gtdbtk.bac120.ar53.summary.tsv',required=True)
    ars('-c',dest="qs",action="store",help='MetaFlowX_all_Original_Bins_all_level_quality.xls',required=True)
    ars('-o',dest="outfile",action="store",help='outfile, default:./',required=True)
    ars('-d',dest="dastool",action="store",help='all sample dastools fa name info',required=True)

    return vars(p.parse_args())

def paste_QS_file(qsF,qsDir):
    with open(qsF,'r') as qsFile:
        header = qsFile.readline()
        for q in  qsFile:
            ql = q.strip().split('\t')
            qsDir[ql[0]] = [ql[1],ql[2],ql[6],ql[8],ql[-1]] 

def paste_gtdb_file(gtdbF,gtdbDir):
    with open(gtdbF,'r') as gtdbFile:
        header = gtdbFile.readline()
        for g in  gtdbFile:
            gl = g.strip().split('\t')
            gtdbDir[gl[0]] = [gl[1],gl[-1]]
            
def get_fastani(summaryF,fastaniDir):
    with open(summaryF,'r') as summaryFile:
        header = summaryFile.readline()
        for s in summaryFile:
            sl = s.strip().split("\t")
            fastaniDir[sl[0]] = sl[2]

def get_bin_sample_info(dastoolFile, bin_sample_Dir):
    with open(dastoolFile,'r') as dastoolF:
        for d in dastoolF:
            dl = d.strip().split("\t")
            binfa = dl[1].rstrip(".fa")
            bin_sample_Dir[binfa] = dl[0]

def merge_table(renameF,outF,qsDir,gtdbDir):
    with open(renameF,'r') as renameFile, open(outF,'w') as outFile:
        outFile.write("BinID\toriginalID\tCompleteness\tContamination\tContig_N50\tGenome_Size\tQS\tGTDB_taxonomy\tNCBI_taxonomy\tfastani_reference\tnative_sample\n")
        for b in renameFile:
            bl = b.strip().split(' ')
            oldBinID = bl[0].rstrip(".fa")
            newBinID = bl[1].rstrip(".fa")
            if oldBinID in qsDir and newBinID in gtdbDir and newBinID in fastaniDir and oldBinID in bin_sample_Dir:
                qs_txt = "\t".join(qsDir[oldBinID])
                taxo_txt = "\t".join(gtdbDir[newBinID])
                outFile.write(f"{newBinID}\t{oldBinID}\t{qs_txt}\t{taxo_txt}\t{fastaniDir[newBinID]}\t{bin_sample_Dir[oldBinID]}\n")
            else:
                print(f"{oldBinID}\t{newBinID}\n")
            
def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,qsDir,gtdbDir,fastaniDir,bin_sample_Dir
    pars = ParsReceiver()
    qsDir,gtdbDir,fastaniDir,bin_sample_Dir = {},{},{},{}
    paste_QS_file(pars['qs'],qsDir)
    paste_gtdb_file(pars['gtdb'],gtdbDir)
    get_fastani(pars['summary'],fastaniDir)
    get_bin_sample_info(pars['dastool'],bin_sample_Dir)
    merge_table(pars['map'],pars['outfile'],qsDir,gtdbDir)
    

if __name__ == '__main__':
    main()