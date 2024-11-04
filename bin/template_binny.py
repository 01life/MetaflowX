#!/usr/bin/env python

import argparse
import sys


text = """mem:
  # If your HPC resource has a high memory capacity node you can set this to
  # the number of high-memory cores you want to use  (set to 0 otherwise)
  # and specify the amount of memory per core (e.g. if a node has 260 gb of
  # RAM and 10 cores it would be 26).
  big_mem_avail: 0
  big_mem_per_core_gb: 26
  # Memory per core of your computing resource.
  normal_mem_per_core_gb: 8
# Path to a temporary directory to write to.
tmp_dir: tmp
raws:
  # Path to test assembly fasta.
  assembly:  "{fa}"
  # Path to a test bam files to calculate depth from.
  metagenomics_alignment: "{bam}"
  # Path to an average depth per contig tsv file. Leave empty for initialization.
  contig_depth: ""
# Sample name
sample: "test"
# Path to desired output dir binny should create and store results in.
outputdir: "{output_path}"
# Absolute path to put binny dbs in. If left empty they will be put into 'database' in the binny main dir.
db_path: ""
# If you want to use existing environments containing the environments for Snakemake, Prokka
# and/or Mantis you can either input the absolute paths to yaml files or env names here.
# Otherwise, leave empty and binny will take care of the installations.
# If you already have Snakemake in your path set: snakemake_env="in_path".
snakemake_env: "base"
prokka_env: "base"
mantis_env: "base"
# Set path for conda envs to be installed to. By default, they will be put in `conda` in the binny dir.
conda_source: ""
# Input a list, e.g. '2,3,4'.
kmers: '2,3,4'
# Mask potentially disruptive contig regions (e.g. rRNA and CRISPR elements) from k-mer counting
mask_disruptive_sequences: 'True'
# Extract single contig MAGs of at least 90% purity and 92.5% completeness
extract_scmags: 'True'
# Will use coassembly mode, starting with contigs >= 500 bp instead of high threshold, decreasing, if set to 'on' or
# if 'auto' and multiple depth files are detected
# Choose between: 'auto', 'on' , 'off''
coassembly_mode: 'auto'
# Binny prefilters assemblies based on N<X> value to try and take as much information as possible into account,
# while minimizing the amount of noise. Be aware that, depending on the assembly quality, low values as the N<X>
# might results in leaving out a large proportion of the assembly (if the max_cont_length cutoffs are set high as well).
NX_value: 90
# Minimum and maximum contig length. Caps value from NX filtering.
min_cont_length_cutoff: 2250  # 500
max_cont_length_cutoff: 2250  # 1000
# Minimum and maximum length of contigs containing CheckM markers. Caps value from NX filtering (initial value is NX_value / 2).
min_cont_length_cutoff_marker: 2250  # 250
max_cont_length_cutoff_marker: 2250  # 400
# Maximum number of contigs binny uses. If the number of available
# contigs after minimum size filtering exceeds this, binny will
# increase the minimum size threshold until the maximum is reached.
# Prevents use of excessive amounts of memory on large assemblies.
# Default should ensure adequate performance, adjust e.g. according
# to available memory.
max_n_contigs: 5.0e5  # 3.5e5
# Maximum marker set lineage depth to check bin quality with:
# 0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'
max_marker_lineage_depth_lvl: 2
# Distance metric for opentSNE and HDBSCAN.
distance_metric: 'manhattan'
embedding:
  # Maximum number of binny iterations.
  max_iterations: 50
clustering:
  # Increasing the HDBSCAN cluster selection epsilon beyond 0.5
  # is not advised as it might massively increase run time, but it might
  # help recover fragmented genomes that would be missed with lower settings.
  hdbscan_epsilon_range: '0.250,0.000'
  # Adapted from the HDBSCAN manual: 'Measure of how conservative the
  # clustering should be. With larger values, more points will be declared
  # as noise, and clusters will be restricted to progressively more dense areas.'.
  hdbscan_min_samples_range: '1,5,10' # '1,2,4,6,8'
  # Use depth as additional dimension during the initial clustering.
  include_depth_initial: 'False'
  # Use depth as additional dimension during the main clusterings.
  include_depth_main: 'False'
bin_quality:
  # Minimum value binny will lower completeness to while running.
  min_completeness: 72.5
  # Completeness threshold binny wilt begin with.
  start_completeness: 92.5
  # Minimum purity for bins to be selected.
  purity: 95
write_contig_data: 'True'
"""


p = argparse.ArgumentParser(description='Output the Binny YAML text to the standard output stream.')
p.add_argument('-f', help='FA file path.')
p.add_argument('-b', help='BAM file path.')
p.add_argument('-o', help='Binny results output path.')
a = p.parse_args()
sys.stdout.write(text.format(fa=a.f, bam=a.b, output_path=a.o))

