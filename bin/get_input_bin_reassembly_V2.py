#!/usr/bin/env python3

import argparse as ap
from collections import defaultdict
import logging
from pathlib import Path
import sys
from typing import List, Dict
import pandas as pd


def check_file_can_be_opened(file: Path):
    """
    Check if a file can be opened.

    Args:
        file (Path): The path to the file to be checked.

    Raises:
        FileNotFoundError: If the file does not exist.
        IsADirectoryError: If the file is not a regular file.
        IOError: If the file is not readable.
    """
    if not file.exists():
        raise FileNotFoundError(f"{file} does not exist.")
    if not file.is_file():
        raise IsADirectoryError(f"{file} is not a file.")
    try:
        with open(file, "r"):
            pass
    except:
        raise IOError(f"{file} is not readable.")
    

def read_bins_rename_map_file(file: Path) -> Dict[str, str]:
    """
    Reads a bins_rename_map file and returns a dictionary containing the information.

    Args:
        file (Path): The path to the bins_rename_map file.

    Returns:
        Dict[str, str]: A dictionary where the keys are the new ID of bins 
        and the values are the raw names of bins.
    """
    bins_rename_map = {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split()
            if len(items) != 2:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            bins_rename_map[items[1].split('.fa')[0]] = items[0].split('.fa')[0]
    return bins_rename_map


def read_fastq_paths_file(file: Path) -> tuple[str, List[Path]]:
    """
    Reads a file containing paths to FASTQ files and returns a dictionary 
    mapping sample IDs to their corresponding FASTQ file paths.
    
    Args:
        file (Path): The path to the file containing the FASTQ paths.
        
    Returns:
        tuple[str, List[Path]]: A tuple containing the sample ID and 
        a list of Path objects representing the FASTQ file paths.
    """
    sample_fastq_paths_map = {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith(("#", "id")):
                continue
            items = line.split(",")
            if len(items) < 3:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            elif len(items) == 3:
                sample_id, fq1, fq2 = items
            else:
                sample_id, fq1, fq2, contig = items
            sample_fastq_paths_map[sample_id] = (Path(fq1), Path(fq2))
    return sample_fastq_paths_map


def read_bins_quality_file(file: Path) -> Dict[str, float]:
    """
    Reads a bins quality file and returns a dictionary mapping bin names to quality values.

    Args:
        file (Path): The path to the bins quality file.

    Returns:
        Dict[str, float]: A dictionary mapping bin names (str) to quality values (float).

    Raises:
        ValueError: If the bins quality file is invalid.
    """
    bins_quality_map = {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith(("#", 'Name')):
                continue
            items = line.split("\t")
            if len(items) != 13:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            bin_name = items[0]
            quality = float(items[-1])
            bins_quality_map[bin_name] = quality
    return bins_quality_map


def read_bins_quality_file_2(file: Path) -> Dict[str, float]:
    """
    Extracts "BinID" as the key and "QS" column as value from a list of dictionaries.

    Args:
        data (list of dictionaries): The data to extract from.

    Returns:
        dict: A dictionary with "BinID" as keys and corresponding "QS" values.
    """

    # Load data into a DataFrame
    df = pd.read_csv(file,sep="\t")

    # Create a dictionary with BinID as key and QS as value
    bins_quality_map = df.set_index('BinID')['QS'].to_dict()

    return bins_quality_map


def read_bins_abundance_file(file: Path, remove_samples: List[str]=[]) -> tuple[Dict[str, List[float]], List[str]]:
    """
    Reads a bins abundance file and returns a dictionary mapping bin IDs to a list of abundances.
    
    Args:
    - file (Path): The path to the bins abundance file.
    - remove_samples (List[str], default=[]): A list of sample names to be removed from the abundance calculations.
        
    Returns:
    - bins_abundance_map (Dict[str, List[float]]): A dictionary mapping bin IDs to a list of abundances.
    - sample_names (List[str]): A list of sample names.
    """
    bins_abundance_map = {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split("\t")
            if len(items) == 1:
                logging.error(f"\n{file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            if line.startswith("BinID"):
                # remove  postfix "_bin_abundance" ".sorted RPKM" ".sorted Mean"
                # sample_names = [sample_name.split("_bin_abundance")[0] for sample_name in items[1:]]
                sample_names = [sample_name.split(".sorted Mean")[0] for sample_name in items[1:]]
                indices = list(range(len(sample_names)))
                if len(remove_samples) > 0:
                    indices = [i for i, sample_name in enumerate(sample_names) if sample_name not in remove_samples]
                    sample_names = [sample_names[i] for i in indices]
                continue
            if len(items) != len(sample_names) + len(remove_samples) + 1:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            bin_id = items[0]
            bin_abundances = [float(item) for i, item in enumerate(items[1:]) if i in indices]
            bins_abundance_map[bin_id] = bin_abundances
    return bins_abundance_map, sample_names


def read_gtdb_summary_file(file: Path) -> tuple[Dict[str, str], Dict[str, str]]:
    """
    Reads a GTDB summary file and returns two dictionaries: gtdb_bins_reference_map and gtdb_species_reference_map.

    Parameters:
        file (Path): The path to the GTDB summary file.

    Returns:
        tuple[Dict[str, str], Dict[str, str]]: A tuple containing the two dictionaries:
        - gtdb_bins_reference_map: A dictionary mapping bin IDs to reference IDs.
        - gtdb_species_reference_map: A dictionary mapping species names to reference IDs.
    """
    gtdb_bins_reference_map, gtdb_species_reference_map = {}, {}
    with open(file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith(("#", "user_genome")):
                continue
            items = line.split("\t")
            if len(items) != 20:
                logging.error(f"\nline {i+1} in {file} is invalid.\n")
                raise ValueError(f"{file} is invalid.")
            bin_id, species, reference_name = items[:3]
            gtdb_bins_reference_map[bin_id] = reference_name
            species = species.split(';')[-1]
            gtdb_species_reference_map[species] = reference_name
    return gtdb_bins_reference_map, gtdb_species_reference_map


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


def get_input_by_bins_quality(
        bin_genomes_dir: Path, 
        # bin_rename_file: Path, 
        fastq_paths_file: Path, 
        bin_quality_file: Path, 
        bin_abundance_file, 
        gtdb_summary_file: Path, 
        gtdb_genome_paths_file: Path,
        quality_threshold: float,
        abundance_threshold: float,
        present_threshold: float,
        remove_samples: List[str]=[]
    ) -> List:
    """
        Retrieves input files and generates a list of selected bins for further reassembly.

        Args:
            bin_genomes_dir (Path): Path to the directory containing bin files.
            # bin_rename_file (Path): Path to the file containing bin rename map.
            fastq_paths_file (Path): Path to the file containing fastq paths map.
            bin_quality_file (Path): Path to the file containing bin quality map.
            bin_abundance_file (Path): Path to the file containing bin abundance map.
            gtdb_summary_file (Path): Path to the file containing GTDB summary.
            gtdb_genome_paths_file (Path): Path to the file containing GTDB genome paths.
            remove_samples (List[str], optional): List of samples to remove. Defaults to [].
        
        Returns:
            List: A list of selected bins for further reassembly.
    """
    bin_genomes = bin_genomes_dir.glob('*/*.fa')
    # bins_rename_map = read_bins_rename_map_file(bin_rename_file)
    fastq_paths_map = read_fastq_paths_file(fastq_paths_file)
    # bins_quality_map = read_bins_quality_file(bin_quality_file)
    bins_quality_map = read_bins_quality_file_2(bin_quality_file)
    if len(remove_samples) > 0:
        bins_abundance_map, sample_names = read_bins_abundance_file(bin_abundance_file, remove_samples=remove_samples)
    else:
        bins_abundance_map, sample_names = read_bins_abundance_file(bin_abundance_file)
    if len(sample_names) == 0:
        logging.error(f"No samples found in {bin_abundance_file} excluding --remove-samples.")
        raise ValueError(f"No samples found in {bin_abundance_file} excluding --remove-samples.")
    gtdb_bins_reference_map, _ = read_gtdb_summary_file(gtdb_summary_file)
    gtdb_genome_paths_map = read_gtdb_genome_paths_file(gtdb_genome_paths_file)
    
    results = []
    for bin_genome in bin_genomes:
        bin_id = bin_genome.stem
        if bin_id not in bins_quality_map:
            logging.warning(f"{bin_id} is not in {bin_quality_file}. Skipping.")
            continue
        quality = bins_quality_map[bin_id]

        if bin_id not in bins_abundance_map:
            logging.warning(f"{bin_id} is not in {bin_abundance_file}. Skipping.")
            continue
        abundances = bins_abundance_map[bin_id]
        if bin_id not in gtdb_bins_reference_map:
            logging.warning(f"{bin_id} is not in {gtdb_summary_file}. Skipping.")
            continue
        reference_id = gtdb_bins_reference_map[bin_id]
        if reference_id not in gtdb_genome_paths_map:
            logging.warning(f"{reference_id} is not in {gtdb_genome_paths_file}. Skipping.")
            continue
        ref_genome = gtdb_genome_paths_map[reference_id]
        if not ref_genome.exists():
            raise FileNotFoundError(f"Reference genome: {ref_genome} does not exist.")

        # Select the sample FASTQ with the highest abundance.
        max_abundance_index = abundances.index(max(abundances))
        max_abundance_sample = sample_names[max_abundance_index]
        if max_abundance_sample not in fastq_paths_map:
            raise FileNotFoundError(f"Fastq of {max_abundance_sample} does not exist.")
        max_abundance_fq1 = fastq_paths_map[max_abundance_sample][0]
        max_abundance_fq2 = fastq_paths_map[max_abundance_sample][1]

        # Filter bins with QS > 90 and present (abundance > 1%) > 10% for further reassembly.
        present = len([i for i in abundances if i > abundance_threshold]) / len(abundances)
        if quality > quality_threshold and present > present_threshold:
            logging.info(f"{bin_id} is selected. QS={quality}, P={present}")
            fq1_paths, fq2_paths = [], []
            match_sample = []
            remove_samples = set(remove_samples)
            for i, sample in enumerate(sample_names):
                # Do not use samples with very low abundance for reassembly.
                if abundances[i] < abundance_threshold:
                    logging.warning(f"{bin_id}: {sample} has low abundance {abundances[i]}. Skipping.")
                    continue
                if sample not in remove_samples:
                    if sample not in fastq_paths_map:
                        raise FileNotFoundError(f"Fastq of {sample} does not exist.")

                    fq1_paths.append(fastq_paths_map[sample][0])
                    fq2_paths.append(fastq_paths_map[sample][1])
                    match_sample.append(sample)
            results.append([bin_id, bin_genome, ref_genome, max_abundance_sample, ','.join(map(str, match_sample))])


    return results


def get_input_target_bins(
        input_file: Path, 
        fastq_paths_file: Path,
        bin_abundance_file: Path,
        gtdb_summary_file: Path, 
        gtdb_genome_paths_file: Path, 
        remove_samples: List[str]=[]
    ) -> List:
    """
    Generate the input target bins based on the given input file, fastq paths file, GTDB summary file, and GTDB genome paths file.

    Args:
        input_file (Path): The path to the input file.
        fastq_paths_file (Path): The path to the fastq paths file.
        gtdb_summary_file (Path): The path to the GTDB summary file.
        gtdb_genome_paths_file (Path): The path to the GTDB genome paths file.
        remove_samples (List[str], optional): A list of samples to be removed. Defaults to [].

    Returns:
        List: A list of lists containing the generated input target bins.
            Each element in the inner list contains the bin ID, bin genome, reference genomes, and fastq paths.
    """
    fastq_paths_map = read_fastq_paths_file(fastq_paths_file)
    if len(remove_samples) > 0:
        bins_abundance_map, sample_names = read_bins_abundance_file(bin_abundance_file, remove_samples=remove_samples)
    else:
        bins_abundance_map, sample_names = read_bins_abundance_file(bin_abundance_file)
    if len(sample_names) == 0:
        logging.error(f"No samples found in {bin_abundance_file} excluding --remove-samples.")
        raise ValueError(f"No samples found in {bin_abundance_file} excluding --remove-samples.")
    _, gtdb_species_reference_map = read_gtdb_summary_file(gtdb_summary_file)
    gtdb_genome_paths_map = read_gtdb_genome_paths_file(gtdb_genome_paths_file)

    results = []
    with open(input_file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split("\t")
            if len(items) != 4:
                logging.error(f"\nline {i+1} in {input_file} is invalid.\n")
                raise ValueError(f"{input_file} is invalid.")
            bin_id, bin_genome, ref_genomes, species = items
            # When specifying reference genomes, multiple can be allowed, separated by commas.
            if ref_genomes != "null":
                ref_genome_paths = [Path(path) for path in ref_genomes.split(",")]
                exist_flags = [path.exists() for path in ref_genome_paths]
                for i, flag in enumerate(exist_flags):
                    if not flag:
                        logging.warning(f"{bin_id}: {ref_genome_paths[i]} does not exist. Skipping.")
                        continue
            #  When specifying species, look up the corresponding reference in the GTDB database.
            elif species != "null":
                if species not in gtdb_species_reference_map:
                    logging.warning(f"{bin_id}: {species} is not in {gtdb_summary_file}. Skipping.")
                    continue
                reference_id = gtdb_species_reference_map[species]
                if reference_id not in gtdb_genome_paths_map:
                    logging.warning(f"{bin_id}: {species} is not in {gtdb_genome_paths_file}. Skipping.")
                    continue
                ref_genomes = gtdb_genome_paths_map[reference_id]
                if not ref_genomes.exists():
                    raise FileNotFoundError(f"Reference genome: {ref_genomes} does not exist.")
            else:
                logging.error(f"\nline {i+1} in {input_file} is invalid. Only one of Species and Reference Genomes can be null\n")
                raise ValueError(f"{input_file} is invalid.")
            
            # Select the sample FASTQ with the highest abundance.
            abundances = bins_abundance_map[bin_id]
            max_abundance_index = abundances.index(max(abundances))
            max_abundance_sample = sample_names[max_abundance_index]
            if max_abundance_sample not in fastq_paths_map:
                raise FileNotFoundError(f"Fastq of {max_abundance_sample} does not exist.")
            max_abundance_fq1 = fastq_paths_map[max_abundance_sample][0]
            max_abundance_fq2 = fastq_paths_map[max_abundance_sample][1]
            # Based on the input FASTQ path file, find FASTQ paths excluding removed samples.
            fq1_paths, fq2_paths = [], []
            match_sample =[]
            remove_samples = set(remove_samples)
            for i, sample in enumerate(sample_names):
                # Do not use samples with very low abundance for reassembly.
                if abundances[i] < abundance_threshold:
                    logging.warning(f"{bin_id}: {sample} has low abundance {abundances[i]}. Skipping.")
                    continue
                if sample not in remove_samples:
                    if sample not in fastq_paths_map:
                        raise FileNotFoundError(f"Fastq of {sample} does not exist.")

                    fq1_paths.append(fastq_paths_map[sample][0])
                    fq2_paths.append(fastq_paths_map[sample][1])
                    match_sample.append(sample)
            # results.append([bin_id, bin_genome, ref_genomes, max_abundance_fq1, max_abundance_fq2, ','.join(map(str, fq1_paths)), ','.join(map(str, fq2_paths))])
            results.append([bin_id, bin_genome, ref_genomes, max_abundance_sample, ','.join(map(str, match_sample))])

    return results


def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-----------------------------------------------------------------------------------------------------------------

    Generate the input file which can be used for preprocess bin reassembly.\n
    Output file format: bin_id\\tbin_genome\\tref_genomes\\tfq1_paths\\tfq2_paths\n\n
    Example1(set the --target_bins_file):\n
    \tpython3 get_input_bin_reassembly.py 
        \t-i target_bins.txt -b bin_dir  
        \t--fq_paths_file read.csv 
        \t--gtdb_summary_file gtdbtk.bac120.summary.tsv 
        \t--gtdb_genome_paths_file genome_paths.tsv 
        \t-o bin_reassembly.txt\n\n
    Example2(only set the --bin_dir):\n
    \tpython3 get_input_bin_reassembly.py 
        \t-b bin_dir 
        \t--fq_paths_file read.csv 
        \t--gtdb_summary_file gtdbtk.bac120.summary.tsv 
        \t--gtdb_genome_paths_file genome_paths.tsv 
        \t-o bin_reassembly.txt

-----------------------------------------------------------------------------------------------------------------'''
)
    ars = p.add_argument
    ars( "-i",
        "--target_bins_file",
        type=str,
        help="Path to the input file of target bins, each line contains:\n\tbin_id\\tbin_genome\\tref_genome2,ref_genome2,...\\tspecies\n"
        "If the ref_genomes set to null, the species will be used to find the reference genome.",)

    ars(
        "-b",
        "--bin_dir",
        type=str,
        help="Direactory to the result of bin refinement.",
        required=True
    )
    ars(
        "--fq_paths_file",
        type=str,
        help="Path to the file containing the fastq paths.",
        required=True,
    )
    ars(
        "--gtdb_summary_file",
        type=str,
        help="Path to the file containing the GTDB summary.",
        required=False,
    )
    ars(
        "--bin_quality_file",
        type=str,
        help="Path to the file containing the bins checkm2 result.",
        required=False,
    )
    ars(
        "--gtdb_genome_paths_file",
        type=str,
        help="Path to the file containing the GTDB genome paths.[$database/release214/fastani/genome_paths.tsv]",
        required=True,
    )
    ars(
        "--remove_samples",
        type=str,
        nargs="*",
        help="A list of samples to be removed.",
        default=[],
    )
    ars(
        "-o",
        "--output_file",
        type=str,
        help="Path to the output file.",
        required=True,
    )
    ars(
        "-q",
        "--quality_score",
        type=float,
        default=90,
        help="Minimum bin quality score to be selected. [%(default)s].",
    )
    ars(
        "-a",
        "--abundance",
        type=float,
        default=0.01,
        help="Minimum bin abundance threshold. [%(default)s].",
    )
    ars(
        "-p",
        "--present",
        type=float,
        default=0.1,
        help="Reassembly bin is only performed if bin_present(bin_abundance > --abundance) > --present. [%(default)s].",
    )
    ars(
        "--prefix",
        type=str,
        default='MetaFlowX',
        help=" The pipeline prefix, must same as the nextflow config setting. '--pipeline_prefix' .default:[MetaFlowX].",
    )

    return vars(p.parse_args())

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,abundance_threshold,bin_dir,pipeline_prefix

    pars = ParsReceiver()
    pipeline_prefix = pars['prefix']

    logging.basicConfig(
        level=logging.DEBUG,
        format=" %(asctime)s - %(levelname)s - %(message)s",
    )
    logging.debug("Start of program.")

    bin_dir = Path(pars['bin_dir']).resolve()
    fastq_path_file = Path(pars['fq_paths_file']).resolve()

    if pars['gtdb_summary_file']:
        gtdb_summary_file = Path(pars['gtdb_summary_file']).resolve()
    else:
        # gtdb_summary_file = bin_dir / "05.BinSet" / "052.Annotation" / "GTDB" / "bin_QS_taxonomy_summary.xls"
        gtdb_summary_file = Path(bin_dir / "05.BinSet" / "052.Annotation" / "GTDB" / "gtdbtk.bac120.ar53.summary.tsv")


    remove_samples = pars['remove_samples']
    output_file = Path(pars['output_file']).resolve()
    gtdb_genome_paths_file = Path(pars['gtdb_genome_paths_file']).resolve()

    quality_score_threshold = pars['quality_score']
    abundance_threshold = pars['abundance'] * 100
    present_threshold = pars['present']

    if not bin_dir.exists():
        logging.error(f"Bin directory {bin_dir} does not exist.")
        raise FileNotFoundError(f"{bin_dir} does not exist.")
    check_file_can_be_opened(fastq_path_file)
    check_file_can_be_opened(gtdb_summary_file)
    check_file_can_be_opened(gtdb_genome_paths_file)
    
    logging.info(f"Bin directory: {bin_dir}")
    logging.info(f"fastq path file: {fastq_path_file}")
    logging.info(f"gtdb summary file: {gtdb_summary_file}")
    logging.info(f"gtdb genome paths file: {gtdb_genome_paths_file}")
    logging.info(f"remove samples: {remove_samples}")
    logging.info(f"output file: {output_file}")
    
    # If the user specifies bins to be reassembled, directly read the input file and extract FASTQ paths and reference genomes based on the information.
    if pars['target_bins_file']:
        target_bins_file = Path(pars['target_bins_file']).resolve()
        # bin_abundance_file = bin_dir / "06.BinsetProfile" / "061.BinAbundance" / "MetaFlowX_bins_AvgDepth.xls"
        bin_abundance_file = Path(bin_dir / "06.BinsetProfile" / "061.BinAbundance" / (pipeline_prefix + "_CoverM_bins_mean_rename.xls"))
        gtdb_summary_file = Path(bin_dir / "05.BinSet" / "052.Annotation" / "GTDB" / "gtdbtk.bac120.ar53.summary.tsv")
            
        check_file_can_be_opened(target_bins_file)
        check_file_can_be_opened(bin_abundance_file)
        logging.info(f"Selected target bins from {target_bins_file}.")
        results = get_input_target_bins(
            target_bins_file, 
            fastq_path_file,
            bin_abundance_file,
            gtdb_summary_file, 
            gtdb_genome_paths_file, 
            remove_samples
        )
    # If the user does not specify bins to be reassembled, select bins from the directory with QS > quality_score_threshold (90) & present (abundance > abundance_threshold (1%)) > present_threshold (10%).
    elif pars['bin_dir']:
        bin_genome_dir = Path( bin_dir / "05.BinSet" / "051.UniqueBin" / "HQUniqueBins" )
        bin_rename_file = Path( bin_dir / "05.BinSet" / "051.UniqueBin" / (pipeline_prefix +"_HQ_unique_bins_rename_map.xls"))
        bin_abundance_file = Path( bin_dir / "06.BinsetProfile" / "061.BinAbundance" / ( pipeline_prefix +"_CoverM_bins_mean_rename.xls"))


        if pars['bin_quality_file']:
            bin_quality_file = Path(pars['bin_quality_file']).resolve()
        else:
            bin_quality_file = Path( bin_dir / "05.BinSet" / "052.Annotation" / "GTDB" / "bin_QS_taxonomy_summary.xls" ) #MetaFlowX_all_Original_Bins_all_level_quality.xls

        check_file_can_be_opened(bin_rename_file)
        check_file_can_be_opened(bin_quality_file)
        check_file_can_be_opened(bin_abundance_file)      
        logging.info(f"Selected bins from {bin_genome_dir}.")
        results = get_input_by_bins_quality(
            bin_genome_dir, 
            # bin_rename_file, 
            fastq_path_file, 
            bin_quality_file, 
            bin_abundance_file,
            gtdb_summary_file,
            gtdb_genome_paths_file,
            quality_score_threshold,
            abundance_threshold,
            present_threshold,
            remove_samples
        )
    else:
        raise ValueError("Either target_bins_file or bin_dir must be specified.")
    
    with open(output_file, 'w') as o:
        # o.write("#BinID\tBinGenome\tRefGenomes\tHighestAbundanceFq1\tHighestAbundanceFq2\tFq1\tFq2\n")
        o.write("#BinID\tBinGenome\tRefGenomes\tHighestAbundance\tFq\n")
        for result in results:
            o.write('\t'.join(map(str, result)) + '\n')
    logging.info(f"End of program.\n\n")



if __name__ == '__main__':
    main()