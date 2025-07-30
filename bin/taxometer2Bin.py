#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import re

def read_taxonomy_mapping(tsv_path: Path) -> dict:
    """Read TSV file and group contigs by taxonomy."""
    tax2contigs = {}
    
    with open(tsv_path,'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            last_taxo = parts[1].split(";")[-1]
            tax2contigs.setdefault(last_taxo, []).append(parts[0])
    return tax2contigs

def split_fasta_by_taxonomy(fasta_path: Path, tax2contigs: dict, prefix: str, outdir: Path):
    """Split fasta by taxonomy group and write to separate bin files."""
    contig2seq = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    n=0

    with open("tomoma_bin_taxonomy_info.tsv",'w') as bin_taxoF:
        bin_taxoF.write(f"BinID\tcontig_taxonomy\n")
        for tax, contigs in tax2contigs.items():
            bin_seqs = []
            for c in contigs:
                if c in contig2seq:
                    bin_seqs.append(contig2seq[c])
                else:
                    print(f"Warning: contig '{c}' not found in FASTA.")

            if bin_seqs:
                outfile = outdir / f"{prefix}_{n}.fa"
                SeqIO.write(bin_seqs, outfile, "fasta")

                bin_taxoF.write(f"{prefix}_{n}\t{tax}\n")
                print(f"✅ Wrote {len(bin_seqs)} contigs to {outfile}")
                n+=1

def main():
    parser = argparse.ArgumentParser(description="Split contig FASTA into taxonomy bins.")
    parser.add_argument("-i", "--taxometer", required=True, help="Taxometer output TSV file")
    parser.add_argument("-f", "--fasta", required=True, help="Contig FASTA file")
    parser.add_argument("-p", "--prefix", default="taxometerSpecies", help="Prefix for output bin files")
    parser.add_argument("-o", "--outdir", default="bins", help="Output directory")
    args = parser.parse_args()

    contig_taxonomy_path = Path(args.taxometer)
    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tax2contigs = read_taxonomy_mapping(contig_taxonomy_path)
    split_fasta_by_taxonomy(fasta_path, tax2contigs, args.prefix, outdir)

if __name__ == "__main__":
    main()
