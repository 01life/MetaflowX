#!/usr/bin/env python

import argparse as ap
from collections import defaultdict
import logging
from pathlib import Path
import sys
from typing import List, Dict
import pandas as pd


def get_union(*lists):
    return list(set().union(*lists))

def get_target_bin_suitable_sample(df: pd.DataFrame, bin_id: str, min_cut: float = 1.0) -> list:
    """
    Extract sample names (column names) suitable for a target bin based on a coverage threshold.

    Args:
        df (pd.DataFrame): Input data frame containing coverage information
        bin_id (str): Target bin ID
        min_cut (float): Minimum coverage threshold. Defaults to 1.0

    Returns:
        list: List of sample names (column names) with coverage values greater than the threshold for the target bin

    Raises:
        ValueError: If the specified bin_id is not in the data frame
    """
    if bin_id not in df.index:
        raise ValueError(f"Specified bin_id '{bin_id}' not found in the data")
    row = df.loc[bin_id]
    suitable_samples = row[row > min_cut].index.tolist()
    return suitable_samples

def read_gtdb_genome_paths_file(file: Path) -> Dict[str, Path]:
    """
    Read a GTDB genome paths file and return a dictionary mapping reference IDs to genome paths.

    Parameters:
        file (Path): The path to the GTDB genome paths file.

    Returns:
        Dict[str, Path]: A dictionary mapping reference IDs to genome paths.
    """
    base_dir = file.parent
    gtdb_genome_paths_map = {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split()
            if len(items) != 2:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            reference_file, genome_path = items
            reference_id = reference_file.split("_genomic.fna")[0]
            ref_genome = base_dir / genome_path / reference_file
            gtdb_genome_paths_map[reference_id] = ref_genome
    return gtdb_genome_paths_map


def get_bin_native_reference(bin_renameFile:str, min_completeness:float=90, min_contamination:float = 5, min_qs:float=65) -> tuple[Dict[str, str], Dict[str, str],List[str]]:
    bin_native_sampe_dir={}
    bin_fastani_reference_dir = {}
    target_bin_list=[]
    with open(bin_renameFile,'r') as bin_renameF:
        header  = bin_renameF.readline()
        for n in bin_renameF:
            binID,originalID,completeness,contamination,contig_N50,genome_Size,qs,gtdb_taxonomy,ncbi_taxonomy,fastani_reference,native_sample  = n.strip().split("\t")
            # sampleID = originalID.strip().split("contigs2bin.tsv_")[1].split(".")[0]
            sampleID = native_sample
            bin_native_sampe_dir[binID] = sampleID
            #if fastani_reference != "N/A":
            bin_fastani_reference_dir[binID] = fastani_reference
            
            if float(completeness) < min_completeness or float(contamination) > min_contamination or float(qs) < min_qs:
                target_bin_list.append(binID)

    return(bin_native_sampe_dir,bin_fastani_reference_dir,target_bin_list)

def get_target_sample(countFile:Path, 
                      meanFile:Path,
                      target_bin_smaple_file:Path,
                      bin_list: List[str]=[], 
                      min_count:float=10000, 
                      min_cov:float=1, 
                      singleAssembly: bool = False)  :
    
    cov = pd.read_csv(meanFile, sep='\t', index_col=0, dtype={0: str})
    count = pd.read_csv(countFile, sep='\t', index_col=0, dtype={0: str})

    with open(target_bin_smaple_file,'w') as outF:
        # Run the function
        for onebin in bin_list:
            one_cov_list = get_target_bin_suitable_sample(cov, onebin, min_cov)
            one_count_list = get_target_bin_suitable_sample(count, onebin, min_count)
            oneSample_list = get_union(one_cov_list,one_count_list)
            Sample_list_txt = ",".join(oneSample_list)

            if bin_fastani_reference_dir[onebin] != 'N/A':
                ref_genome = gtdb_genome_paths_Dir[bin_fastani_reference_dir[onebin]]
            else:
                ref_genome = 'N/A'

            if singleAssembly:
                outF.write(f'{onebin}\t{bin_native_sampe_dir[onebin]}\t{ref_genome}\t{bin_native_sampe_dir[onebin]}\n')
            else:
                outF.write(f'{onebin}\t{Sample_list_txt}\t{ref_genome}\t{bin_native_sampe_dir[onebin]}\n')



def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-----------------------------------------------------------------------------------------------------------------

    select bin ID and target sample ID to reassembly.\n
    Output file format: bin_id\\tbin_genome\\tref_genomes\\tfq1_paths\\tfq2_paths\n\n
    Example1(co-assembly):\n
    \tpython3 get_input_bin_reassembly.py 
        \t-i bin_QS_taxonomy_summary.xls
        \t-c MetaFlowX_CoverM_bins_count_rename.xls
        \t-m MetaFlowX_CoverM_bins_trimmed_mean_rename.xls\n\n
    Example2(signle assembly):\n
    \tpython3 get_input_bin_reassembly.py 
        \t-i bin_QS_taxonomy_summary.xls
        \t-c MetaFlowX_CoverM_bins_count_rename.xls
        \t-m MetaFlowX_CoverM_bins_trimmed_mean_rename.xls\n

-----------------------------------------------------------------------------------------------------------------'''
)
    ars = p.add_argument
    ars( "-i",
        "--bin_summary",
        type=str,
        help="Path to bin_QS_taxonomy_summary.xls",
        required=True
    )

    ars(
        "-c",
        "--count",
        type=str,
        help="Path to MetaFlowX_CoverM_bins_count_rename.xls",
        required=True
    )
    ars(
        "-m",
        "--mean",
        type=str,
        help="Path to MetaFlowX_CoverM_bins_trimmed_mean_rename.xls",
        required=True
    )
    ars(
        "--minCompleteness",
        type=float,
        default=90,
        help="Minimum bin Completeness threshold. [%(default)s].",
    )
    ars(
        "--minContamination",
        type=float,
        default=5,
        help="Minimum bin Contamination threshold. [%(default)s].",
    )
    ars(
        "--minQS",
        type=float,
        default=65,
        help="Minimum bin quality score to be selected. [%(default)s].",
    )

    ars(
        "--minCount",
        type=float,
        default=10000,
        help="Minimum bin reads count threshold. [%(default)s].",
    )
    ars(
        "--minDepth",
        type=float,
        default=1,
        help="Minimum bin depth threshold. [%(default)s].",
    )
    ars(
        "--singleAssembly",
        action='store_true',
        help="Use co-assembly pattern by default. Set '--singleAssembly' to use only native sample map reads for assembly. Default is False.",)

    ars(
        "--prefix",
        type=str,
        default='MetaFlowX',
        help=" The pipeline prefix, must same as the nextflow config setting. '--pipeline_prefix' .default:[MetaFlowX].",
    )
    ars(
        "--gtdb_genome_paths_file",
        type=str,
        help="Path to the file containing the GTDB genome paths.[$database/release214/fastani/genome_paths.tsv]",
        required=True,
    )

    return vars(p.parse_args())

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,bin_native_sampe_dir,bin_fastani_reference_dir,gtdb_genome_paths_Dir

    pars = ParsReceiver()
    pipeline_prefix = pars['prefix']

    gtdb_genome_paths_file = Path(pars['gtdb_genome_paths_file']).resolve()

    gtdb_genome_paths_Dir = read_gtdb_genome_paths_file(gtdb_genome_paths_file)


    bin_native_sampe_dir,bin_fastani_reference_dir,target_bin_list  = get_bin_native_reference(pars['bin_summary'],pars['minCompleteness'],pars['minContamination'],pars['minQS'])

    get_target_sample(pars['count'],pars['mean'],f'{pipeline_prefix}_pick4optimize_bin_sample.txt',target_bin_list,pars['minCount'],pars['minDepth'],pars['singleAssembly'])



if __name__ == '__main__':
    main()