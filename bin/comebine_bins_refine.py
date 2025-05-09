#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
from collections import Counter

def parse_args():
	parser = argparse.ArgumentParser(description="Refine bins from multiple bin folders. Author:lianglifeng V1	2023-09-27")
	parser.add_argument('-c', '--contig', required=True, help="Paths to contigs fa.")
	parser.add_argument('-info', required=True, help="Paths to contigs squeence information file which including squences length.")
	parser.add_argument('-i','--sample', required=True, help="sampleID")
	parser.add_argument('-bins', nargs='+', required=True, help="Paths to binner contigs2bin.tsv file.")
	parser.add_argument('-o', '--output', required=False,default=os.getcwd(), help="Output directory.")
	parser.add_argument('-ms', '--min_size', type=float, default=0.5, help="Minimum size of refined bins in Mbp (default=0.5 Mbp).")
	parser.add_argument('-mx', '--max_size', type=float, default=20, help="Maximum size of refined bins in Mbp (default=20 Mbp).")
	parser.add_argument('-t', '--threads', type=int, default=os.cpu_count(), help="Number of threads to use (default: all available)")
	return parser.parse_args()

def create_output_dir(output_dir):
	os.makedirs(output_dir, exist_ok=True)

def parse_fats(contigF, contigDir):
	for record in SeqIO.parse(contigF, "fasta"):
		contigDir[record.id] =  record.seq


def get_contig_lengths(contig_info_file):
	len_dir = {}
	
	with open(contig_info_file,'r') as cf:
		header = cf.readline()
		for contig in cf:
			cl = contig.strip().split("\t")
			len_dir[cl[0]] =  int(float(cl[1]))
	return(len_dir)


def process_contigs2bin_file(contig2binFile,binner_name):
	
	contig_bin_dir = {}
	binner_bin_contig_dir = {}

	with open(contig2binFile,'r') as cbF:
		for b in cbF:
			contig_id, bin_id = b.strip().split("\t")

			contig_bin_dir.setdefault(binner_name,{})[contig_id] = bin_id

			binner_bin_contig_dir.setdefault(binner_name, {}).setdefault(bin_id, []).append(contig_id)

	
	return(contig_bin_dir, binner_bin_contig_dir)


def refine_bins(binner_combination, total_contig_bin_dir, total_binner_bin_contig_dir, sorted_contig_lengths, min_size, max_size, sampleID):
	binning_contig_list = set()  # Use a set to avoid duplicates.
	combine_bin = {}
	multi_binner_common_bin = {}
	binner_list = binner_combination
	
	n = 0
	for contig in sorted_contig_lengths:
		if contig not in binning_contig_list:
			bin_contig_list = []

			# Directly process multiple binners to reduce duplicate merging of lists.
			for binner in binner_list:

				binid = total_contig_bin_dir.get(binner, {}).get(contig)

				if binid is None:
					print(f"Warning: {contig} not found for binner {binner}")
				else:
					contig_match_bins_contig_list = total_binner_bin_contig_dir[binner][binid]
					bin_contig_list.extend(contig_match_bins_contig_list)
			
			# Count occurrences and extract common contigs.
			counter = Counter(bin_contig_list)
			common_contig_list = [item for item, count in counter.items() if count == len(binner_list)]

			binning_contig_list.update(common_contig_list)

			# Only add to combine_bin when common_contig_list is not empty.
			if common_contig_list:
				refined_bin_size = sum(sorted_contig_lengths[contig] for contig in common_contig_list) / 1e6  # Convert to Mbp
				if refined_bin_size >= min_size and refined_bin_size <= max_size:
					combine_bin[n] = common_contig_list
					n += 1
	newbinid = str(sampleID)+"_"+":".join(binner_combination)
	multi_binner_common_bin[newbinid] = combine_bin

	return multi_binner_common_bin


def process_bin(binner_combination, binid, common_bin_contig_dir, output_dir, contigDir,sample):
	"""Processes a single bin by writing contigs and their sequences to files.

	Args:
		binner_combination (str): Binner combination.
		binid (str): Bin ID.
		common_bin_contig_dir (dict): Common bin contig directory.
		output_dir (str): Output directory.
		contigDir (dict): Dictionary containing contig sequences.
	"""

	# Generate bin-specific identifiers and file paths
	binID = f"{binner_combination}.{binid}"

	mulit_bin = binner_combination.split("_")[-1].split(".")[0]
	create_output_dir(os.path.join(output_dir,mulit_bin))
	onebinFa_path = os.path.join(output_dir,mulit_bin, f"{binID}.fa")

	# Open FASTA file for writing contigs


	with open(onebinFa_path, 'w') as onebinFa, open(os.path.join(output_dir, f"{sample}_allcontigs2bin.txt"), 'a') as contigbinF:
		for contig in common_bin_contig_dir[binid]:
			# Write to .contigs2bin.tsv
			filename = binner_combination.replace(":", "__")
			with open(os.path.join(output_dir, f"{filename}.contigs2bin.tsv"), 'a') as outfile:
				outfile.write(f"{contig}\t{binID}\n")
			
			# Write contig sequence to FASTA file
			onebinFa.write(f">{contig}\n{contigDir[contig]}\n")
			
			#save contig:[binid ]
			contigbinF.write(f"{contig}\t{binID}\n")

			

	return f"{binID} processing complete."

def write_refined_bins(refined_bins, output_dir, contigDir, threads, sample):
	"""Writes refined bins to files, incorporating multi-processing with ProcessPoolExecutor.

	Args:
		refined_bins (dict): Dictionary containing refined bins.
		output_dir (str): Output directory.
		contigDir (dict): Dictionary containing contig sequences.
		threads (int): Number of threads for multi-processing.
	"""

	# Ensure output directory exists
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	# Create a list to hold futures
	futures = []

	# Create a pool of workers with ProcessPoolExecutor
	with ProcessPoolExecutor(max_workers=threads) as executor:
		for binner_combination in refined_bins:
			common_bin_contig_dir = refined_bins[binner_combination]

			# Submit jobs to the executor
			for binid in common_bin_contig_dir:
				futures.append(executor.submit(
					process_bin, binner_combination, binid, common_bin_contig_dir, output_dir, contigDir,sample))

		# Wait for all jobs to complete and retrieve their results
		for future in as_completed(futures):
			print(future.result())


def process_single_biner_output(binner_list,total_binner_bin_contig_dir,output_dir,contigDir,sample):
	n = 0
	for binner_name, bin_contig in total_binner_bin_contig_dir.items():
		for oldBinid, one_bin_contig in bin_contig.items():

			binID = f"{sample}_{binner_name}:onlyOneBiner.{n}"
			
			n += 1
			create_output_dir(os.path.join(output_dir,f"{binner_name}:onlyOneBiner"))
			onebinFa_path = os.path.join(output_dir,f"{binner_name}:onlyOneBiner", f"{binID}.fa")

			# Open FASTA file for writing contigs
			with open(onebinFa_path, 'w') as onebinFa, open(os.path.join(output_dir, f"{sample}_allcontigs2bin.txt"), 'a') as contigbinF, open(os.path.join(output_dir, f"{sample}_{binner_name}__onlyOneBiner.contigs2bin.tsv"), 'a') as outfile:
				for contig in one_bin_contig:
					# Write to .contigs2bin.tsv
					outfile.write(f"{contig}\t{binID}\n")
					
					# Write contig sequence to FASTA file
					onebinFa.write(f">{contig}\n{contigDir[contig]}\n")
					
					#save contig:[binid ]
					contigbinF.write(f"{contig}\t{binID}\n")
					

def main():
	args = parse_args()
	if len(args.bins) < 2:
		print("Only get one Binner result, MetaflowX will not combine bins.")

		# raise ValueError("At least two bin folders are required.")
	
	create_output_dir(args.output)

	contigDir = {}
	parse_fats(args.contig, contigDir)

	# 1 # get all contig len
	contig_lengths = get_contig_lengths(args.info)
	sorted_contig_lengths = dict(sorted(contig_lengths.items(), key=lambda item: item[1], reverse=True))


	# 2 # get contig and bins relationship dictionary
	binner_list  = []
	total_contig_bin_dir  = {}
	total_binner_bin_contig_dir = {}
	total_combine_refine_bin = {}

	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		futures = []
		for contig2bin in args.bins:
			binner_name = os.path.basename(contig2bin).strip().split(".contigs2bin")[0]
			binner_list.append(binner_name)
			print(binner_name)
			futures.append(executor.submit(process_contigs2bin_file, contig2bin, binner_name))
		
		for future in as_completed(futures):
			contig_bin_dir, binner_bin_contig_dir = future.result()
			total_contig_bin_dir.update(contig_bin_dir)

			total_binner_bin_contig_dir.update(binner_bin_contig_dir)

	if len(binner_list) > 1:

		# 3 # get binner combination and merge bin
		all_binner_combinations = []
		for r in range(2, len(binner_list) + 1):
			combinations = list(itertools.combinations(binner_list, r))
			all_binner_combinations.extend(combinations)

		
		print(all_binner_combinations)
		with ProcessPoolExecutor(max_workers=args.threads) as executor:
			futures = []
			for bin_combination in all_binner_combinations:
				futures.append(executor.submit(refine_bins, bin_combination, total_contig_bin_dir, total_binner_bin_contig_dir, sorted_contig_lengths, args.min_size, args.max_size, args.sample))

			for future in as_completed(futures):
				combine_refine_bin = future.result()
				total_combine_refine_bin.update(combine_refine_bin)

		write_refined_bins(total_combine_refine_bin, args.output,contigDir,args.threads,args.sample)

	else:
		print(binner_list)
		process_single_biner_output(binner_list,total_binner_bin_contig_dir,args.output,contigDir,args.sample)
		# print(contig_bin_dir)
		# print(total_contig_bin_dir)
		# print(total_binner_bin_contig_dir)


if __name__ == "__main__":
	main()