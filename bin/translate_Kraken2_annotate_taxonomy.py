#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path

def parse_taxonomy_file(taxo_file):
    df = pd.read_csv(taxo_file, sep='\t')
    taxodir = {}
    for line in df['TaxonomyLevel']:
        if isinstance(line, str) and 's_' in line:
            species = [x for x in line.split('|') if x.startswith('s_')]
            if species:
                taxodir[species[0].replace('s_', '').strip()] = line
    return taxodir

def parse_kraken_output(kraken_file, taxodir):
    result = []
    with open(kraken_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            contig = fields[1]
            species_info = fields[2].split()
            if len(species_info) < 2:
                continue
            species_name = f"{species_info[0]} {species_info[1]}"
            taxonomy = taxodir.get(species_name, '')
            result.append((contig, taxonomy))
    return result

def main():
    parser = argparse.ArgumentParser(description="Map Kraken2 output to full taxonomy level.")
    parser.add_argument('--taxonomy', required=True, help="Bracken species taxonomy file (e.g. sample0_bracken_species_mpa.xls)")
    parser.add_argument('--kraken', required=True, help="Kraken2 output file (e.g. sample0_output.txt)")
    parser.add_argument('--output', default='taxonomy.tsv', help="Output file name (default: taxonomy.tsv)")
    args = parser.parse_args()

    taxodir = parse_taxonomy_file(args.taxonomy)
    mapped = parse_kraken_output(args.kraken, taxodir)

    with open(args.output, 'w') as out:
        out.write("contigs\tpredictions\n")
        for contig, taxonomy in mapped:
            if taxonomy != "":
                out.write(f"{contig}\t{taxonomy.replace('|', ';')}\n")

if __name__ == "__main__":
    main()
