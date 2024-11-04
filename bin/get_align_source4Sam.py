#!/usr/bin/env python

import sys
from collections import defaultdict as dfd


def multi_read(fns):
    for fn in fns:
        for line in open(fn):
            yield line


def find_best_align_and_output(fns, names, opt):
    names = set(names)
    ans = dfd(lambda: (0, '-'))
    it = multi_read(fns) if type(fns) is list else fns
    for line in it:
        if line[0] == '@': continue
        sps = line.rstrip('\n').split('\t')
        name, flag, source = sps[:3]
        if names and name not in names: continue  # Restricted, skip.
        if source != '*':
            if int(flag) & 0x100: # Mapped to multiple locations.
                score = int(sps[4])
                if score > ans[name][0] and score != 255:
                    ans[name] = (score, source)
            else:  # Mapped to only one location.
                opt.write(name+'\t'+source+'\n')
        # else:  # Not mapped.
        #     output.write(name+'\t'+'-'+'\n')
    for i in ans:
        opt.write(i+'\t'+ans[i][1]+'\n')


if __name__ == '__main__':
    if '-h' in sys.argv[1:]: exit('Find the highest quality matching gene names in the SAM file.\n[read1, read2, ...]\tRead names, defaults to all.\n-f\tFile of read names separated by newlines.\n-s\tPath to the SAM file (can input multiple). Defaults to standard input stream.\n-o\tOutput path. Defaults to standard output stream.')
    input_ = sys.stdin
    output = sys.stdout
    args = {j: [] for j in ['reads', '-f', '-s', '-o']}
    status = 'reads'
    for arg in sys.argv[1:]:
        if arg in args:
            status = arg
        else:
            args[status].append(arg)
    if args['-f']:
        for line in open(args['-f'][0]):
            args['reads'].append(line.rstrip('\n'))
    if args['-s']:
        input_ = args['-s']
    if args['-o']:
        output = open(args['-o'][0], 'w')
    args['reads'] = [i.encode() for i in args['reads']]
    find_best_align_and_output(input_, args['reads'], output)
