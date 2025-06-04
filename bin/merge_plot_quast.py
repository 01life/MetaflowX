#!/usr/bin/env python

import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Merge QUAST report TSV files and generate interactive Plotly visualizations.")
parser.add_argument('--result_dir', type=str, required=True, help="Directory containing *_report.tsv files")
parser.add_argument('--prefix', type=str, required=False, default='MetaflowX', help="Prefix for output file")
args = parser.parse_args()

# Set working directory from command-line argument
working_dir = args.result_dir
merged_tsv = os.path.join(working_dir, f"{args.prefix}_quast_report.txt")
output_html = os.path.join(working_dir, f"{args.prefix}_quast_report.html")

# Validate working directory
if not os.path.isdir(working_dir):
    print(f"Error: {working_dir} is not a valid directory.")
    exit(1)

# Find all *_report.tsv files in the working directory
tsv_files = [f for f in os.listdir(working_dir) if f.endswith('_report.tsv')]

# Initialize an empty list to store DataFrames
dfs = []

# Read and process each TSV file
for file in tsv_files:
    # Read TSV file, setting first column as index
    df = pd.read_csv(os.path.join(working_dir, file), sep='\t', index_col=0)
    # Extract sample name from filename (remove '_report.tsv')
    sample_name = file.replace('_report.tsv', '')
    # Rename column to sample name
    df.columns = [sample_name]
    dfs.append(df)

# Merge all DataFrames column-wise
if dfs:
    merged_df = pd.concat(dfs, axis=1)
    
    # Save merged DataFrame to TSV in the working directory
    merged_df.to_csv(merged_tsv, sep='\t')
    
    # Transpose DataFrame for plotting
    merged_df_T = merged_df.T
    
    # Select core metrics for plotting
    metrics = ['N50', 'Total length', '# contigs', 'Largest contig', 'GC (%)']
    # Ensure metrics exist in the DataFrame
    available_metrics = [m for m in metrics if m in merged_df_T.columns]
    
    # Create subplots for each metric
    fig = make_subplots(
        rows=len(available_metrics),
        cols=1,
        subplot_titles=available_metrics,
        vertical_spacing=0.1
    )
    
    # Plot each metric
    for i, metric in enumerate(available_metrics, 1):
        # Convert metric values to numeric, coercing errors to NaN
        y_values = pd.to_numeric(merged_df_T[metric], errors='coerce')
        # Create bar plot for each metric
        fig.add_trace(
            go.Bar(
                x=merged_df_T.index,  # Sample names
                y=y_values,  # Numeric metric values across samples
                name=metric
            ),
            row=i,
            col=1
        )
        
        # Update y-axis title
        fig.update_yaxes(title_text=metric, row=i, col=1)
    
    # Update layout
    fig.update_layout(
        height=300 * len(available_metrics),
        width=800,
        title_text="QUAST Assembly Metrics Comparison",
        showlegend=False
    )
    
    # Save plot as HTML in the working directory
    fig.write_html(output_html)
    
    print(f"Merged TSV saved as '{merged_tsv}'")
    print(f"Interactive plot saved as '{output_html}'")
else:
    print(f"No *_report.tsv files found in {working_dir}.")