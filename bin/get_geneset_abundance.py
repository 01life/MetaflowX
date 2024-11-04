#!/usr/bin/env python

# -*- encoding: utf-8 -*-

import re, argparse
import os, gzip
from multiprocessing import Process, Pipe
from threading import RLock, Thread


class PSet:
    """
    Process Set
    A shared collection between processes, connected using the add_link function.
    add_link returns a closure.
    """

    def __init__(self):
        self.pipes = []
        self.s = set()
        self.rlock = RLock()
        self.finish_signal = -1

    def start(self, pipe_n):

        def th(pipe_n):
            while True:
                rc = pipe_n.recv()
                if rc == self.finish_signal: break
                if rc in self.s:
                    pipe_n.send(True)
                else:
                    pipe_n.send(False)
                    self.rlock.acquire()
                    self.s.add(rc)
                    self.rlock.release()

        Thread(target=th, args=(pipe_n,), daemon=True).start()

    def add_link(self):
        p1, p2 = Pipe()
        self.pipes.append(p2)
        self.start(p2)

        return p1


def read_gene_length(fn):
    #  rt: 22 s, rb: 16 s
    f = open(fn, 'rb')
    fl = f.readline().rstrip(b'\n').lower().split(b'\t')
    index_id, index_name, index_length = fl.index(b'id'), fl.index(b'name'), fl.index(b'length')
    # print('debug id, name. length:', index_id, index_name, index_length)
    marker2glen = {}
    for line in f:
        sps = line.rstrip(b'\n').split(b'\t')
        marker2glen[sps[index_name]] = [sps[index_id], float(sps[index_length])]
    return marker2glen


def multi_read(fns, mode='rb', start_point=0):
    fsizes = [os.path.getsize(i) for i in fns]
    # Calculate the start file and start index.
    for i in fsizes:
        if start_point > i:
            start_point -= i
            fns = fns[1:]
            continue
        else:
            break

    def get():
        for index, fn in enumerate(fns):
            f = open(fn, mode=mode)
            if not index: f.seek(start_point)
            for line in f:
                yield line

    return get


def stat_sam(fns, m2gl, query_pipe, ans_pipe, file_start_index=0, next_file_start_index=0, insert_size=600,
             se_count=1.0, pe_count=0.5):
    get_line = multi_read(fns, 'rb', file_start_index)
    ans = {}

    first = True
    for line in get_line():
        if file_start_index > next_file_start_index:  # The first newline after reading to the next point.
            break
        if first and file_start_index:
            first = False
            file_start_index += len(line)
            continue
        else:
            first = False
        file_start_index += len(line)

        if line[0] == 64:  # ord('@') == 64
            continue
        sps = line.split(b'\t')  # splits

        marker = sps[2].split(b'/')[0]
        read = sps[0].split(b'/')[0]  # Remove the trailing “/1/2”.

        if marker not in m2gl:
            # print('debug: unkonwn marker:', marker)
            continue  # Skip unknown markers.

        if sps[2].endswith(b'*'):
            # print('debug: unmapped reads')
            continue  # Skip reads that did not map successfully.

        count = 0.0
        flag = int(sps[1])

        if flag & 0x100:
            # print('debug: multi_mapping reads：flag=', bin(flag))
            continue  # Skip reads that mapped to multiple locations.

        is_pe_mapping = flag & 0x1  # PE
        if not is_pe_mapping:
            count = se_count
        else:
            this_reverse = flag & 0x10  # This read matches the negative strand.
            paired_reverse = flag & 0x20  # The other read matches the negative strand.

            if flag & 0x2 and (
                    (this_reverse and not paired_reverse) or (not this_reverse and paired_reverse)):  # situation 3
                count = pe_count
            else:
                _id, glen = m2gl[marker]
                base_index = int(sps[3])
                # CIGAR(column 6) -> read_len
                read_len = sum([int(i) for i in re.findall(rb'\d+', sps[5])])

                is_in_right_end = glen - base_index < insert_size  # read in left side of reference strand
                is_in_left_end = base_index < insert_size - read_len  # right side of ref

                if flag & 0x8:  #The other read cannot be mapped.
                    if (this_reverse and is_in_left_end) or (not this_reverse and is_in_right_end):  # situation 1
                        # only aligned in one strand
                        count = se_count
                else:
                    # 5' ->; 3' <-
                    if (not this_reverse and is_in_right_end) or (this_reverse and is_in_left_end):
                        if sps[6] == b'=':
                            query_pipe.send(read)
                            inside = query_pipe.recv()
                            if not inside:  # situation 5; only count one time for each
                                count = se_count
                        else:  # situation 6
                            count = se_count
        ans[marker] = ans.get(marker, 0.0) + count
    query_pipe.send(-1)
    ans_pipe.send(ans)
    return ans


def single_thread_stat_sam(fns, m2gl, insert_size=600, se_count=1.0, pe_count=0.5):
    get_line = multi_read(fns, 'rb')
    reads_set = set()
    ans = {}

    for line in get_line():
        if line[0] == 64:  # ord('@') == 64
            continue
        sps = line.split(b'\t')  # splits

        marker = sps[2].split(b'/')[0]
        read = sps[0].split(b'/')[0]  # Remove the trailing “/1/2”.

        if marker not in m2gl:
            # print('debug: unkonwn marker:', marker)
            continue  # Skip unknown markers.

        if sps[2].endswith(b'*'):
            # print('debug: unmapped reads')
            continue  # Skip reads that did not map successfully.

        count = 0.0
        flag = int(sps[1])

        if flag & 0x100:
            # print('debug: multi_mapping reads：flag=', bin(flag))
            continue  # Skip reads that mapped to multiple locations.

        is_pe_mapping = flag & 0x1  # PE
        if not is_pe_mapping:
            count = se_count
        else:
            this_reverse = flag & 0x10  # This read matches the negative strand.
            paired_reverse = flag & 0x20  # The other read matches the negative strand.

            if flag & 0x2 and ((this_reverse and not paired_reverse) or (not this_reverse and paired_reverse)):  # situation 3
                count = pe_count
            else:
                _id, glen = m2gl[marker]
                base_index = int(sps[3])
                # CIGAR(column 6) -> read_len
                read_len = sum([int(i) for i in re.findall(rb'\d+', sps[5])])

                is_in_right_end = glen - base_index < insert_size  # read in left side of reference strand
                is_in_left_end = base_index < insert_size - read_len  # right side of ref

                if flag & 0x8:   # The other read cannot be mapped.
                    if (this_reverse and is_in_left_end) or (not this_reverse and is_in_right_end):  # situation 1
                        # only aligned in one strand
                        count = se_count
                else:
                    # 5' ->; 3' <-
                    if (not this_reverse and is_in_right_end) or (this_reverse and is_in_left_end):
                        if sps[6] == b'=':
                            if read not in reads_set:  # situation 5; only count one time for each
                                count = se_count
                                reads_set.add(read)
                        else:  # situation 6
                            count = se_count

        ans[marker] = ans.get(marker, 0.0) + count
    return ans


def compute_abundance(
        sam_filenames,
        database_fn,
        sample_name=None,
        output_folder=None,
        prefix=None,
        insert_size=600,
        se_count=1.0,
        pe_count=0.5,
        use_gzip=True,
        process_num=1,
        show_zero_abundance=True,
):
    # print('debug: reading database')
    m2gl = read_gene_length(database_fn)

    if process_num > 1:
        total_fsize = sum([os.path.getsize(i) for i in sam_filenames])
        mission_splits = [(int(total_fsize / process_num * i), int(total_fsize / process_num * (i + 1))) for i in
                          range(process_num)]
        # print('debug: mission splits:', mission_splits)

        pset = PSet()
        pros = []
        pipes = []
        for index, idx in enumerate(mission_splits):
            query_pipe = pset.add_link()
            ans_pipe, p2 = Pipe()
            pipes.append(p2)
            pros.append(
                Process(
                    target=stat_sam,
                    args=(
                        sam_filenames,
                        m2gl,
                        query_pipe,
                        ans_pipe,
                        idx[0],
                        idx[1],
                        insert_size,
                        se_count,
                        pe_count,
                    ),
                    daemon=True,
                )
            )

        for i in pros:
            i.start()

        ans = {}
        sub_ans = []
        for index, i in enumerate(pipes):
            sub_ans.append(i.recv())

        keys = []
        for i in sub_ans:
            keys += list(i.keys())
        keys = set(keys)
        for i in keys:
            ans[i] = sum([j.get(i, 0.0) for j in sub_ans])
    else:  # single thread
        ans = single_thread_stat_sam(
            sam_filenames,
            m2gl,
            insert_size=insert_size,
            se_count=se_count,
            pe_count=pe_count,
        )

    # mkdir
    output_folder = output_folder or './'
    if not os.path.exists(output_folder): os.mkdir(output_folder)
    prefix = prefix or 'sample'
    output_filename_ab = os.path.join(output_folder, prefix + '_abundance' + ('.xls.gz' if use_gzip else '.xls'))
    output_filename_total_ab = os.path.join(output_folder, prefix + '_total_abundance.xls')

    output = (gzip.open if use_gzip else open)(output_filename_ab, 'wb')

    # write header
    output.write(b'\t' + sample_name.encode() + b'\n')
    # compute abundance and output
    total_ab = 0.0
    id2ab = {}
    for marker in m2gl:
        reads_num = ans.get(marker, 0.0)
        if reads_num == 0.0 and not show_zero_abundance: continue
        _id, glen = m2gl[marker]
        ab = reads_num / glen
        id2ab[_id] = ab
        total_ab += ab

    for _id in sorted(id2ab, key=int):
        output.write(_id + b'\t' + str(id2ab[_id] and (id2ab[_id] / total_ab)).encode() + b'\n')

    # output total abundance
    open(output_filename_total_ab, 'w').write(f'ID\ttotal_abundance\n{sample_name}\t{total_ab}\n')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', nargs='+', required=True, help='Input file names, which could be multiple.')
    p.add_argument('-s', '--sample', default='sample', help='Sample name.')
    p.add_argument('-g', '--gene', required=True, help='gene length infor txt,include 3 columns,1st:ID index[0 strat], 2nd: name, 3rd: length ')
    p.add_argument('-o', '--outpath', default='./', help='Output path.')
    p.add_argument('-p', '--prefix', default='sample', help='The prefix of output file names.')
    p.add_argument('-n', '--nprocess', default=1, type=int, help='Process num.')
    p.add_argument('-sc', '--se-count', default=1.0, type=float, help='Single end count, default to be 1.0.')
    p.add_argument('-pc', '--pe-count', default=0.5, type=float, help='Paired end count, default to be 0.5.')
    p.add_argument('-isize', '--insert-size', dest='isize', default=500, type=int, help='Max insert size.')
    p.add_argument('-nz', '--no-zero', dest='nz', action='store_true', help='Do not show zero abundance results in output.')
    p.add_argument('-ng', '--nogzip', action='store_true', help='Do not use gzip for .pr output file.')
    return vars(p.parse_args())


def main():
    args = parse_args()
    compute_abundance(
        args['input'],
        args['gene'],
        sample_name=args['sample'],
        output_folder=args['outpath'],
        prefix=args['prefix'],
        insert_size=args['isize'],
        se_count=args['se_count'],
        pe_count=args['pe_count'],
        use_gzip=not args['nogzip'],
        process_num=args['nprocess'],
        show_zero_abundance=not args['nz'],
    )


if __name__ == '__main__':
    main()

