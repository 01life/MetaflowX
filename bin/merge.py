#!/usr/bin/env python

import argparse
import os
import pandas as pd
import subprocess, sys
from itertools import takewhile,chain

def merge(sample, files, path, prefix):
    """
    Outputs the table join of the given pre-split string collection.
    :param  sampe:   One or more split lines from which data are read.
    :param  files:       files need to merge
    :param  path:       Output stream to which matched rows are written.

    """
    
    #Output file name for MPA abundance merge.
    taxid = os.path.join(path, prefix+"_abundance_table_with_taxid.xls")
    notaxid = os.path.join(path, prefix+"_abundance_table.xls")

    #Output file name for IGC/Humann abundance merge. 
    f = files[0]
    suffix = os.path.splitext(f)[-1]    
    if(suffix=='.gz'):
        f = f.replace(suffix,"")
    names = os.path.splitext(os.path.basename(f))[0].split("_")[1:]
    names = "_".join(names)
    fname = os.path.join(path, prefix+"_"+names+".xls")

    #Output stream.
    out = None
    #Output file name.
    outFile = None

    listmpaVersion = set()
    profiles_list = []
    merged_tables = pd.DataFrame

    #Iterate through all merged files.
    for f in files:
        headers = []
        compression = None
        suffix = os.path.splitext(f)[-1]
        if(suffix!='.gz'):
            headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), open(f))]
        else:
            compression = "gzip"
        h_length = len(headers)

        #Process MPA abundance files.
        if(h_length>1):
            skip = h_length
            listmpaVersion.add(headers[0])
            names = headers[-1].split('#')[1].strip().split('\t')
            usecols=[0,1,2]
            index_col=[0,1]
            cols = os.path.splitext(os.path.basename(f))[0].replace('_profile', '')
            out = open(taxid, 'w')
            outFile = taxid
        #Process IGC/Humann abundance files.
        else:
            skip = 1
            names = ["", "relative_abundance"]
            usecols = [0,1]
            index_col = 0
            cols = tuple(sample[files.index(f)])
            out = open(fname, 'w')
            outFile = fname

        if len(listmpaVersion)>1:
            print('merge_metaphlan_tables: profiles from differrent versions of MetaPhlAn, please profile your '
                  'samples using the same MetaPhlAn version.\n')
            return
        
        iIn = pd.read_csv(f, sep='\t', skiprows=skip, names=names, usecols=usecols, index_col=index_col, compression=compression)
        profiles_list.append(pd.Series(data=iIn["relative_abundance"], index=iIn.index, name=cols))
        
    merged_tables = pd.concat([pd.concat(profiles_list, axis=1).fillna(0)], axis=1, sort=True).fillna(0)

    if(len(listmpaVersion)!=0):
        out.write(list(listmpaVersion)[0]+'\n')
        
    merged_tables.to_csv(out, sep='\t')

    if(len(listmpaVersion)!=0):
        out = open(notaxid, 'w')
        newIndex = pd.Index([i[0] for i in merged_tables.index], name=('clade_name'))
        info = pd.DataFrame(merged_tables.values, columns=merged_tables.columns,index=newIndex)
        info.to_csv(out, sep="\t")

    # postProcess(out, outFile)


def postProcess(out, file):

    #Compress merged files larger than 100MB.
    out.seek(0, os.SEEK_END)
    fsize = out.tell()
    f_mb = fsize/float(1024*1024)
    #For testing, temporarily compress files larger than 100MB.
    if(f_mb >= 100):
        try:
            return_code=subprocess.check_call(f'gzip {file}',shell=True)
        except (EnvironmentError, subprocess.CalledProcessError) as e:
            message="Unable to gzip: "+" ".join(f'gzip {file}')
            if hasattr(e, 'output') and e.output:
                message+="\nError message returned from gizp:\n" + e.output
            sys.exit(message)

argp = argparse.ArgumentParser(prog="merge.py",
                               description="Performs a table join on one or more abundance output files.")
argp.add_argument("-i", help="abundance list contails all samples")
argp.add_argument("-p", help="output file prefix", required=False, default='Test')
argp.add_argument("-o", help="path of output file")

def main( ):

    args = argp.parse_args()

    if args.i:
        content = pd.read_csv(open(args.i), delim_whitespace=True, dtype={"id": str})
        ids = [x.strip().split() for x in content["id"]]
        files = [x.strip().split() for x in content["file"]]
        input = list(chain.from_iterable(files))
        merge(ids, input, args.o if args.o else os.getcwd(), args.p)
    else:
        print("missing input argument!")

if __name__ == '__main__':
    main()
