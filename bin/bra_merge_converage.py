#!/usr/bin/env python

import os
import sys
import pandas as pd

def process_files(file_paths):
    data = {}
    sample_names = []
    tad_counter = 0

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            # Skip the first line
            next(file)
            
            tad_counter += 1
            base_name = os.path.basename(file_path).replace(".txt", ".bam")
            sample_names.append(base_name)
            
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) < 5:
                    continue
                
                cn, cl, tad, bam, var = columns
                if cn not in data:
                    data[cn] = {'contigInfo': {'cl': 0, 'tad': 0}, 'bamInfo': {}}
                
                data[cn]['contigInfo']['cl'] = cl
                data[cn]['contigInfo']['tad'] += float(tad)
                data[cn]['bamInfo'][base_name] = {'bam': bam, 'var': var}
    
    return data, sample_names, tad_counter

def print_output(data, sample_names, tad_counter):
    header = ["contigName", "contigLen", "totalAvgDepth"] + \
             [f"{name}\t{name}-var" for name in sample_names]
    print("\t".join(header))
    # sorted_cn = sorted(data.keys(), key=lambda x: (int(x.split('_')[0].replace('k', '')), int(x.split('_')[1])))
    
    for cn in data.keys():
        contig_info = data[cn]['contigInfo']
        bam_info = data[cn]['bamInfo']
        row = [cn, contig_info['cl'], contig_info['tad'] / tad_counter] + \
              [f"{bam_info.get(name, {}).get('bam', '')}\t{bam_info.get(name, {}).get('var', '')}" for name in sample_names]
        print("\t".join(map(str, row)))

if __name__ == "__main__":
    file_paths = sys.argv[1:]
    data, sample_names, tad_counter = process_files(file_paths)
    print_output(data, sample_names, tad_counter)
