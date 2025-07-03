#!/usr/bin/env python

import pandas as pd
import argparse
import os
from typing import List, Dict, Tuple
from pathlib import Path

def parse_taxonomy(tax_string: str) -> Tuple[List[str], List[str]]:
    """Parse taxonomy string into individual levels and cumulative strings.
    
    Args:
        tax_string: Input taxonomy string (e.g., 'd_Bacteria;p_Proteobacteria')
    
    Returns:
        Tuple of (cumulative taxonomy levels, individual taxonomy parts)
    """
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    prefix_to_level = {f"{level[0]}_": level for level in levels}
    tax_parts = tax_string.split(';') if tax_string else []
    tax_levels = []
    individual_parts = []
    current_tax = []

    for part in tax_parts:
        prefix = part[:2] if '_' in part else ''
        if prefix in prefix_to_level:
            level = prefix_to_level[prefix]
            cleaned_part = part[len(prefix):].strip('[]')
            current_tax.append(f"{prefix}{cleaned_part}")
            tax_levels.append(';'.join(current_tax))
            individual_parts.append(cleaned_part)
        else:
            print(f"Warning: Skipping unexpected part '{part}' with prefix '{prefix}'")
    
    return tax_levels, individual_parts

def parse_scores(score_string: str) -> List[float]:
    """Parse scores string into a list of floats.
    
    Args:
        score_string: String of semicolon-separated scores
    
    Returns:
        List of float scores
    """
    return [float(score) for score in score_string.split(';')] if score_string else []

def process_taxonomy_file(input_file: str, output_dir: str, score_threshold: float = 0.95, prefix: str = "All") -> None:
    """Process taxonomy TSV and save separate TSVs for each taxonomic level and deepest levels.
    
    Args:
        input_file: Path to input TSV file
        output_dir: Directory for output TSV files
        score_threshold: Minimum score for inclusion
        prefix: Prefix for output filenames
    """
    # Define taxonomic levels and output files
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    output_dir = Path(output_dir)
    output_files = {level: output_dir / f"{prefix}_{level}.tsv" for level in levels}
    deepest_level_file = output_dir / f"{prefix}_deepest_level.tsv"
    
    # Initialize data storage
    level_data = {level: [] for level in levels}
    deepest_level_data = []
    
    # Read input file
    df = pd.read_csv(input_file, sep='\t', usecols=['contigs', 'predictions', 'scores'])
    
    # Process each row
    for row in df.itertuples(index=False):
        contig, tax_string, score_string = row.contigs, row.predictions, row.scores
        scores = parse_scores(score_string)
        tax_levels, _ = parse_taxonomy(tax_string)
        
        # Find deepest level meeting threshold
        deepest_level_idx = -1
        for i, score in enumerate(scores):
            if score >= score_threshold:
                deepest_level_idx = i
            else:
                break
        
        # Store data for each level up to deepest level
        for i, (level, tax_level) in enumerate(zip(levels, tax_levels)):
            if i <= deepest_level_idx:
                level_scores = ';'.join(map(str, scores[:i+1]))
                level_data[level].append({
                    'contigs': contig,
                    'predictions': tax_level,
                    'scores': level_scores
                })
        
        # Store deepest level data with only contigs and predictions
        if deepest_level_idx >= 0:
            deepest_level_data.append({
                'contigs': contig,
                'predictions': tax_levels[deepest_level_idx]
            })
    
    # Save output files
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save level-specific files
    for level in levels:
        if level_data[level]:
            pd.DataFrame(level_data[level]).to_csv(
                output_files[level], sep='\t', index=False
            )
            print(f"Saved {level}s to {output_files[level]}")
        else:
            print(f"No data for {level}s (all scores below threshold)")
    
    # Save deepest level file with only contigs and predictions
    if deepest_level_data:
        pd.DataFrame(deepest_level_data, columns=['contigs', 'predictions']).to_csv(
            deepest_level_file, sep='\t', index=False
        )
        print(f"Saved deepest levels to {deepest_level_file}")
    else:
        print("No deepest level data (all scores below threshold)")

def main():
    parser = argparse.ArgumentParser(
        description='Process taxonomy TSV file and save separate TSV files per taxonomic level and deepest levels.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input_file', required=True, help='Input TSV file path')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for TSV files')
    parser.add_argument('-p', '--prefix', default="All", help='Output files prefix')
    parser.add_argument('--score_threshold', type=float, default=0.95, help='Score threshold for filtering')

    args = parser.parse_args()
    process_taxonomy_file(args.input_file, args.output_dir, args.score_threshold, args.prefix)

if __name__ == '__main__':
    main()