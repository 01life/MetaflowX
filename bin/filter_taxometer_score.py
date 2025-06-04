#!/usr/bin/env python

import pandas as pd
import argparse
import os

def parse_taxonomy(tax_string):
    """Parse the taxonomy string into individual levels and cumulative strings."""
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    prefix_to_level = {f"{level[0]}_": level for level in levels}  # Map prefixes to levels
    tax_parts = tax_string.split(';') if tax_string else []
    tax_levels = []
    current_tax = []

    for part in tax_parts:
        # Extract the prefix (e.g., 'd_', 'p_')
        prefix = part.split('_')[0] + '_' if '_' in part else ''
        # Check if the prefix is valid
        if prefix in prefix_to_level:
            level = prefix_to_level[prefix]
            # Remove prefix or brackets from the part
            cleaned_part = part[len(prefix):] if part.startswith(prefix) else part.strip('[]')
            current_tax.append(f"{prefix}{cleaned_part}")
            tax_levels.append(';'.join(current_tax))
        else:
            # Handle unexpected prefixes (e.g., second 'd_') by skipping or logging
            print(f"Warning: Skipping unexpected part '{part}' with prefix '{prefix}'")
            continue
    
    return tax_levels
    

def parse_scores(score_string):
    """Parse the scores string into a list of floats."""
    if score_string:
        return [float(score) for score in score_string.split(';')]
    return []

def process_taxonomy_file(input_file, output_dir, score_threshold=0.95, prefix="All"):
    """Process the taxonomy TSV file and save separate TSV files for each taxonomic level."""
    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Define taxonomic levels and corresponding output files
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    output_files = {level: os.path.join(output_dir, f"{prefix}_{level}.tsv") for level in levels}
    
    # Initialize dictionaries to store data for each level
    level_data = {level: [] for level in levels}
    
    # Process each row
    for _, row in df.iterrows():
        contig = row['contigs']
        tax_string = row['predictions']
        scores = parse_scores(row['scores'])
        tax_levels = parse_taxonomy(tax_string)
        
        # Add data for each taxonomic level if score meets threshold
        for i, (level, tax_level) in enumerate(zip(levels, tax_levels)):
            if i < len(scores) and scores[i] >= score_threshold:
                # Collect contig, taxonomy, and scores up to this level
                level_scores = ';'.join(map(str, scores[:i+1]))
                level_data[level].append({'contigs': contig, 'predictions': tax_level, 'scores': level_scores})
    
    # Save each level to its respective TSV file
    os.makedirs(output_dir, exist_ok=True)
    for level in levels:
        if level_data[level]:  # Only create file if there is data
            output_df = pd.DataFrame(level_data[level], columns=['contigs', 'predictions', 'scores'])
            output_df.to_csv(output_files[level], sep='\t', index=False)
            print(f"Saved {level}s to {output_files[level]}")
        else:
            print(f"No data for {level}s (all scores below threshold)")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process taxonomy TSV file and save to separate TSV files per taxonomic level.')
    parser.add_argument('-i', '--input_file', required=True, help='Input TSV file path')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for TSV files')
    parser.add_argument('-p', '--prefix', default="All", help='Output files prefix')
    parser.add_argument('--score_threshold', type=float, default=0.95, help='Score threshold for filtering (default: 0.95)')

    args = parser.parse_args()
    
    # Process the file
    process_taxonomy_file(args.input_file, args.output_dir, args.score_threshold,args.prefix)

if __name__ == '__main__':
    main()