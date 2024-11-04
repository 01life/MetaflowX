#!/usr/bin/env python

import argparse
import os


def recover(fn_tsv, fn_fa, path, binPrefix, _id_idx=0, source_idx=1):
    print(fn_fa, '->', path)
    if not os.path.exists(path): os.mkdir(path)
    fns = set()
    data = {}
    for line in open(fn_tsv):
        sps = line[:-1].split('\t')
        _id, source = sps[_id_idx], sps[source_idx]
        data['>'+_id+'\n'] = source
        fns.add(source)
    s2f = {i: open(os.path.join(path,binPrefix+"."+ i+'.fa'), 'w', encoding='utf-8') for i in fns}
    fr = open(fn_fa)
    ne = next(fr)
    while ne:
        fa = ''
        for line in fr:
            if line and line[0] != '>':
                fa += line
            else:
                if ne in data:
                    s2f[data[ne]].write(ne+fa)
                ne = line
                break
        else:
            if ne in data:
                s2f[data[ne]].write(ne + fa)
            ne = ''


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-t', help='metabinner_res/metabinner_result.tsv ')
    p.add_argument('-f', help='contig fa')
    p.add_argument('-o', help='outpath')
    p.add_argument('-i', help='contig sequence ID start index, default = 0 ', default=0, type=int,required=False)
    p.add_argument('-s', help='bin name start index, default=1', default=1, type=int,required=False)
    p.add_argument('-p', help='bin name start index, default=bin', default="bin",required=False)
    a = p.parse_args()
    recover(a.t, a.f, a.o, a.p, a.i, a.s, )


if __name__ == '__main__':
    main()
