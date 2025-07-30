#!/usr/bin/env python

from pathlib import Path
import argparse
import csv

def load_mapping(mapping_file: Path) -> dict[str, str]:
    """Load contig ID mapping from a two-column file (old ID → new ID)."""
    mapping = {}
    with mapping_file.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                old_id, new_id = parts[:2]
                mapping[old_id] = new_id
    return mapping

def replace_contig_ids(
    input_tsv: Path,
    mapping: dict[str, str],
    output_tsv: Path
) -> None:
    """
    Replace the first column (contig IDs) in the input TSV file
    based on the provided mapping and write the result to a new file.
    """
    with input_tsv.open("r", encoding="utf-8") as infile, \
         output_tsv.open("w", encoding="utf-8", newline="") as outfile:

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        header = next(reader)
        writer.writerow(header)  # Write the header line as is

        for row in reader:
            contig_id = row[0]
            new_id = mapping.get(contig_id, contig_id)  # Keep original if not mapped
            row[0] = new_id
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description="Replace contig IDs in a TSV file using a mapping file"
    )
    parser.add_argument("--map", required=True, type=Path,
                        help="Mapping file (2 columns: old ID, new ID)")
    parser.add_argument("--input", required=True, type=Path,
                        help="Input TSV file with contig IDs to be replaced")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output TSV file with updated contig IDs")
    args = parser.parse_args()

    mapping = load_mapping(args.map)
    replace_contig_ids(args.input, mapping, args.output)

if __name__ == "__main__":
    main()
