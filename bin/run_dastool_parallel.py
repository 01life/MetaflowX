#!/usr/bin/env python

import os
import sys
import glob
import argparse
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args():
    parser = argparse.ArgumentParser(description="Multi-threaded DAS_Tool runner")
    parser.add_argument('--input', required=True, help='Parent directory containing per-sample subfolders')
    parser.add_argument('--cpus', type=int, default=2, help='CPUs per DAS_Tool run (default: 2)')
    parser.add_argument('--dastool_opts', default=' --write_bin_evals --write_bins ', help='Additional options passed to DAS_Tool')
    return parser.parse_args()

def find_sample_dirs(input_dir):
    if os.path.isdir(input_dir):
        return(os.path.abspath(input_dir))

def collect_sample_files(sample_dir):
    sample_id = os.path.basename(sample_dir)
    pattern = os.path.join(sample_dir, f"{sample_id}_*.tsv")
    bin_files = glob.glob(pattern)

    # Find contigs and protein FASTA files in the current directory (./)
    contig_path = next((f for f in glob.glob("./*contigs.fa")), None)
    pep_path = next((f for f in glob.glob("./*protein.fa")), None)

    # Check if files were found
    if contig_path is None:
        raise FileNotFoundError("No contigs file (*.contigs.fa) found in the current directory")
    if pep_path is None:
        raise FileNotFoundError("No protein file (*.protein.fa) found in the current directory")

    binner_map = {}
    for path in bin_files:
        fname = os.path.basename(path)
        if not fname.startswith(f"{sample_id}_") or not fname.endswith(".tsv"):
            continue
        binner = fname.removeprefix(f"{sample_id}_").removesuffix(".tsv")
        binner_map[binner] = os.path.abspath(path)

    if not os.path.exists(contig_path) or not os.path.exists(pep_path):
        print(f"[WARNING] Missing contig or protein for {sample_id}, skipped.")
        return None

    return sample_id, binner_map, contig_path, pep_path


def run_dastool(sample_id, combo, bin_paths, contigs, pep, cpus, dastool_opts):
    combo_label = "_".join(combo)
    label = f"{sample_id}_{combo_label}"
    bin_list = ",".join(bin_paths)

    cmd = [
        "DAS_Tool",
        "-i", bin_list,
        "-c", contigs,
        "-o", label,
        "-p", pep,
        "-t", str(cpus)
    ]
    if dastool_opts:
        cmd += dastool_opts.strip().split()
    
    print(' '.join(cmd))

    try:
        print(f"[INFO] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        # 检查是否生成 bin 结果
        bin_out_dir = f"{label}_DASTool_bins"
        if os.path.isdir(bin_out_dir) and glob.glob(os.path.join(bin_out_dir, "*.fa")):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            fasta2bin_script = os.path.join(script_dir, "Fasta_to_Contig2Bin.sh")
            contig2bin_cmd = f"{fasta2bin_script} -i {bin_out_dir} -e fa > {sample_id}_{combo_label}_contigs2bin.tsv"
            eval_cmd = f"sed 's/^/{sample_id}\\t{combo_label}\\t/g' {label}_allBins.eval > {sample_id}_{combo_label}_allBins.eval.txt"
            subprocess.run(contig2bin_cmd, shell=True)
            subprocess.run(eval_cmd, shell=True)
        else:
            print(f"[INFO] No output bins found for {label}, skipped downstream processing.")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] DAS_Tool failed for {label}: {e}")

def process_sample(sample_dir, cpus, dastool_opts):
    result = collect_sample_files(sample_dir)
    if not result:
        return

    sample_id, binner_map, contigs, pep = result

    if len(binner_map) < 2:
        print(f"[INFO] Less than 2 binners for {sample_id}, skip.")
        return

    combos = []
    for r in range(2, len(binner_map) + 1):
        combos.extend(itertools.combinations(sorted(binner_map.keys()), r))


    for combo in combos:
        bin_paths = [binner_map[b] for b in combo]
        run_dastool(sample_id, combo, bin_paths, contigs, pep, cpus, dastool_opts)
            
def main():
    args = parse_args()
    sample_dir = os.path.abspath(args.input)
    process_sample(sample_dir, args.cpus, args.dastool_opts)


if __name__ == "__main__":
    main()
