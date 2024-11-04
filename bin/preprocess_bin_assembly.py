#!/usr/bin/env python

import argparse
from collections import defaultdict
import logging
import os
import re
from pathlib import Path
import subprocess
import sys
import time
from typing import List
import shutil


def check_file_can_be_opened(file: Path) -> None:
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
    return None


def check_exe_can_be_run(exe: Path) -> None:
    """
    Check if the given executable can be run.

    Args:
        exe (Path): The path to the executable.

    Raises:
        FileNotFoundError: If the executable does not exist.
        PermissionError: If the executable is not executable.
    """
    if not exe.exists():
        raise FileNotFoundError(f"{exe} does not exist.")
    if not os.access(exe, os.X_OK):
        raise PermissionError(f"{exe} is not executable.")
    return None


def get_fq(fq_path):
    path_dir = {}
    for root, _, files in os.walk(fq_path):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith('fq.gz'):
                fileName = os.path.splitext(file)[0]
                sampleID = re.match(r"^(.*)_clean_", fileName).group(1)
                path_dir.setdefault(sampleID,[]).append(Path(file_path).resolve())
    return(path_dir)


def run_cmd(cmd: str, description: str) -> str:
    if description:
        logging.info(f"{description}\n")
    logging.info(f"Running command:\n{cmd}")
    start_time = time.time()
    try:
        result = subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}):\n{}".format(e.cmd, e.returncode, e.stderr.decode()))
    end_time = time.time()
    logging.info(f"Process completed. Time elapsed: {(end_time - start_time) / 60:.5f} miniutes.\n")
    return result.stdout.decode()


def cat_fasta_files(fasta_files: List[Path], output_file: Path) -> None:
    outdir = output_file.parent
    cat_cmd = "cat "
    gunzip_files = []
    for fasta_file in fasta_files:
        if fasta_file.suffix == ".gz":
            gunzip_file = outdir / fasta_file.with_suffix('').name
            gunzip_files.append(gunzip_file)
            gunzip_cmd = f"gzip -d -c {fasta_file} > {gunzip_file}"
            run_cmd(gunzip_cmd, description=f"Decompress {fasta_file}.")
            cat_cmd += f"{gunzip_file} "
        else:
            cat_cmd += f"{fasta_file} "
    cat_cmd += f"> {output_file}"
    run_cmd(cat_cmd, description="Merge reference genomes and bin contigs.")

    # clean up
    for file in gunzip_files:
        file.unlink()
    return None


def cat_fastq_files(fq1_files: List[Path], fq2_files: List[Path], prefix: str) -> tuple[Path, Path]:
    if fq1_files[0].suffix == ".gz":
        fq1_output = Path(prefix + "_1.fq.gz")
        fq2_output = Path(prefix + "_2.fq.gz")
    else:
        fq1_output = Path(prefix + "_1.fq")
        fq2_output = Path(prefix + "_2.fq")

    # merge fastq1
    fq1_cmd = "cat "
    for fastq_file in fq1_files:
        fq1_cmd += f"{fastq_file} "
    fq1_cmd += f"> {fq1_output}"
    # merge fastq2
    fq2_cmd = "cat "
    for fastq_file in fq2_files:
        fq2_cmd += f"{fastq_file} "
    fq2_cmd += f"> {fq2_output}"

    cmd = fq1_cmd + "\n" + fq2_cmd
    run_cmd(cmd, description=f"Merge fastq files of {prefix}.")
    return fq1_output, fq2_output


def bwa_mem_mapping(
    bin_id: str,
    fq1_files: List[Path],
    fq2_files: List[Path],
    ref_fasta: Path,
    outdir: Path,
    bwa_path: str,
    samtools_path: str,
    bwa_mem_args: str,
    samtools_view_args: str,
    threads: int
) -> tuple[List[Path], List[Path]]:
    """
    Executes BWA MEM mapping for multiple FASTQ files to a reference genome.
    
    Args:
        bin_id (str): The ID of the bin.
        fq1_files (List[Path]): A list of Path objects representing the paths to the first paired-end FASTQ files.
        fq2_files (List[Path]): A list of Path objects representing the paths to the second paired-end FASTQ files.
        ref_fasta (Path): A Path object representing the path to the reference genome FASTA file.
        outdir (Path): A Path object representing the output directory where the mapping results will be stored.
        bwa_path (str): The path to the BWA executable.
        samtools_path (str): The path to the Samtools executable.
        bwa_mem_args (str): Additional command-line arguments to be passed to BWA MEM.
        samtools_view_args (str): Additional command-line arguments to be passed to Samtools view.
        threads: The number of threads to use.
        
    Returns:
        tuple[List[Path], List[Path]]: A tuple containing two lists of Path objects representing the paths to the
            generated mapping FASTQ files for each input FASTQ file.
    """
    # Establish BWA index for the new reference genome.
    index_cmd = f"{bwa_path} index {ref_fasta}"
    run_cmd(index_cmd, description="BWA index for reference genome.")

    mapping_fq1_files, mapping_fq2_files = [], []
    for i, fq1_file in enumerate(fq1_files):
        fq2_file = fq2_files[i]
        sample_name = f"sample.{i}"
        bam_file = outdir / f"{bin_id}_{sample_name}.bam"
        mapping_fq1 = outdir / f"{bin_id}_{sample_name}_mapping_1.fq.gz"
        mapping_fq1_files.append(mapping_fq1)
        mapping_fq2 = outdir / f"{bin_id}_{sample_name}_mapping_2.fq.gz"
        mapping_fq2_files.append(mapping_fq2)
        # Execute BWA MEM to map sample FASTQs to the reference genome and use samtools view to filter target reads.
        map_cmd = (f"{bwa_path} mem -t {threads} {bwa_mem_args} {ref_fasta} {fq1_file} {fq2_file} |"
            f"{samtools_path} view -@ {threads} {samtools_view_args} -o {bam_file} -")
        run_cmd(map_cmd, description=f"Mapping reads to reference genome of {bin_id}_{sample_name}.")
        # Execute samtools collate to convert BAM files to FASTQ.
        fastq_cmd = (f"{samtools_path} collate -@ {threads} -u -O {bam_file} |"
                f"{samtools_path} fastq -@ {threads} -1 {mapping_fq1} -2 {mapping_fq2} -N -0 /dev/null -s /dev/null")
        run_cmd(fastq_cmd, description=f"Convert bam to fastq of {bin_id}_{sample_name}.")
    return mapping_fq1_files, mapping_fq2_files


def filter_fq_by_mash(
    bin_id: str,
    mapping_fq1_files: List[Path],
    mapping_fq2_files: List[Path],
    outdir: Path,
    max_abundance_sample_index: int,
    max_dist_threshold: float,
    mash_path: str,
    mash_sketch_args: str,
    threads: int,
) -> tuple[List[Path], List[Path]]:
    """
    Filter fastq files by performing MASH sketch precomputation on the fastq file with the highest abundance and calculating the distance between this file and other fastq files. Returns the filtered fastq files.

    Args:
        bin_id (str): The ID of the bin.
        mapping_fq1_files (List[Path]): A list of Path objects representing the mapping fastq1 files.
        mapping_fq2_files (List[Path]): A list of Path objects representing the mapping fastq2 files.
        outdir (Path): The output directory for the MASH sketch files.
        max_abundance_sample_index (int): The index of the fastq file with the highest abundance in the mapping_fq1_files and mapping_fq2_files lists.
        max_dist_threshold (float): The maximum distance threshold for filtering the fastq files.
        mash_path (str): The path to the MASH executable.
        mash_sketch_args (str): Additional arguments to be passed to the MASH sketch command.
        threads (int): The number of threads to use.

    Returns:
        tuple[List[Path], List[Path]]: A tuple containing two lists of Path objects representing the filtered fastq files. The first list contains the mapping fastq1 files, and the second list contains the mapping fastq2 files.
    """
    if len(mapping_fq1_files) == 1:
        return [mapping_fq1_files,], [mapping_fq2_files,]
    
    # Execute MASH sketch to precompute the FASTQ of the sample with the highest abundance.
    max_abundance_fq1 = mapping_fq1_files.pop(max_abundance_sample_index)
    max_abundance_fq2 = mapping_fq2_files.pop(max_abundance_sample_index)
    if os.path.getsize(max_abundance_fq1) < 1024:
        raise ValueError(f"{max_abundance_fq1} is too small.")
    max_abundance_fq1_msh = outdir / max_abundance_fq1.with_suffix(max_abundance_fq1.suffix + ".msh").name
    max_abundance_fq2_msh = outdir / max_abundance_fq2.with_suffix(max_abundance_fq2.suffix + ".msh").name
    cmd = (f"{mash_path} sketch -o {max_abundance_fq1_msh} {mash_sketch_args} {max_abundance_fq1}\n"
        f" {mash_path} sketch -o {max_abundance_fq2_msh} {mash_sketch_args} {max_abundance_fq2}")
    run_cmd(cmd, description=f"MASH sketching {max_abundance_fq1} and {max_abundance_fq2} which are the highest abundance.")

    dist_file = outdir / f"{bin_id}_dist.txt"
    if dist_file.exists():
        dist_file.unlink()
    for i, mapping_fq1_file in enumerate(mapping_fq1_files):
        # TODO: deprecated
        if os.path.getsize(mapping_fq1_file) < 1024:
            logging.warning(f"{mapping_fq1_file} is too small for reassembly, skipping.")
            continue
        mapping_fq2_file = mapping_fq2_files[i]
        # Execute MASH sketch to precompute the FASTQ of the remaining samples, and use MASH dist to calculate the distance between these FASTQs and the one with the highest abundance.
        mapping_fq1_msh = outdir / mapping_fq1_file.with_suffix(mapping_fq1_file.suffix + ".msh").name
        mapping_fq2_msh = outdir / mapping_fq2_file.with_suffix(mapping_fq2_file.suffix + ".msh").name
        mash_sketch_cmd = (f"{mash_path} sketch -p {threads} -o {mapping_fq1_msh} {mash_sketch_args} {mapping_fq1_file}\n"
                    f" {mash_path} sketch -p {threads} -o {mapping_fq2_msh} {mash_sketch_args} {mapping_fq2_file}")
        run_cmd(mash_sketch_cmd, description=f"MASH sketching {mapping_fq1_file} and {mapping_fq2_file}.")
        mash_dist_cmd = (f"{mash_path} dist -p {threads} {max_abundance_fq1_msh} {mapping_fq1_msh} >> {dist_file}\n"
                    f"{mash_path} dist -p {threads} {max_abundance_fq2_msh} {mapping_fq2_msh} >> {dist_file}")
        run_cmd(mash_dist_cmd, description=f"MASH distance measure between {max_abundance_fq1} and {mapping_fq1_file}.")
    
    # Read the dist file, identify samples with dist < max_dist_threshold, and return the corresponding FASTQs.
    output_fq1_files, output_fq2_files = [max_abundance_fq1, ], [max_abundance_fq2, ]
    i = 0
    fq2_flag = False
    with open(dist_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items = line.split()
            if not fq2_flag:
                fq1_dist = float(items[2])
                # fq1 and fq2 are alternating, indicating that the next line contains the distance results between fq2.
                fq2_flag = True
            else:
                fq2_dist = float(items[2])
                fq2_flag = False
                if fq1_dist < max_dist_threshold and fq2_dist < max_dist_threshold:
                    output_fq1_files.append(mapping_fq1_files[i])
                    output_fq2_files.append(mapping_fq2_files[i])
                i += 1
    return output_fq1_files, output_fq2_files


def preprocess_sample_fastqs(
    bin_id: str,
    fq1_files: List[Path],
    fq2_files: List[Path],
    max_abundance_fq1: Path,
    max_abundance_fq2: Path,
    bin_fasta: Path,
    ref_fastas: List[Path],
    outdir: Path,
    use_single_sample: bool,
    not_use_reference: bool,
    bwa_path: str,
    bwa_mem_args: str,
    samtools_path: str,
    samtools_view_args: str,
    mash_path: str,
    max_dist_threshold: float,
    mash_sketch_args: str,
    threads: int,
) -> tuple[Path, Path]:
    """
    Preprocesses the sample fastq files for a given bin.

    Args:
        bin_id (str): The ID of the bin.
        fq1_files (List[Path]): The list of Path objects representing the first paired-end fastq files.
        fq2_files (List[Path]): The list of Path objects representing the second paired-end fastq files.
        max_abundance_fq1 (Path): The Path object representing the first paired-end fastq file with maximum abundance.
        max_abundance_fq2 (Path): The Path object representing the second paired-end fastq file with maximum abundance.
        bin_fasta (Path): The Path object representing the bin fasta file.
        ref_fastas (List[Path]): The list of Path objects representing the reference fasta files.
        outdir (Path): The Path object representing the output directory.
        bwa_path (str): The path to the BWA executable.
        bwa_mem_args (str): The additional arguments to pass to the BWA-MEM command.
        samtools_path (str): The path to the Samtools executable.
        samtools_view_args (str): The additional arguments to pass to the Samtools view command.
        mash_path (str): The path to the Mash executable.
        max_dist_threshold (float): The maximum distance threshold for filtering samples based on Mash distance.
        mash_sketch_args (str): The additional arguments to pass to the Mash sketch command.

    Returns:
        tuple[Path, Path]: A tuple containing the Path objects representing the output fastq files.
    """
    # If use_single_sample is True, only use the FASTQ of the sample with the highest abundance.
    if use_single_sample:
        fq1_files = [max_abundance_fq1, ]
        fq2_files = [max_abundance_fq2, ]
    else:
        fq1_files.append(max_abundance_fq1)
        fq2_files.append(max_abundance_fq2)
    max_abundance_index = len(fq1_files) - 1
    
    print(
        "\n#############################BWA mapping#####################################\n",
        file=sys.stderr,
    )
    # Extract FASTQs that can be mapped to the reference genome and bins.
    mapping_fq_dir = outdir / "MappingFq"
    mapping_fq_dir.mkdir(exist_ok=True)
    merge_ref_fasta = mapping_fq_dir / f"{bin_id}_ref.fasta"
    if not_use_reference:
        run_cmd(f"ln -sf {bin_fasta} {merge_ref_fasta}")
    else:
        cat_fasta_files(ref_fastas + [bin_fasta,], merge_ref_fasta)
    mapping_fq1_files, mapping_fq2_files = bwa_mem_mapping(
        bin_id,
        fq1_files,
        fq2_files,
        merge_ref_fasta,
        mapping_fq_dir,
        bwa_path,
        samtools_path,
        bwa_mem_args,
        samtools_view_args,
        threads
    )
    print(
        "\n#############################BWA mapping#####################################\n",
        file=sys.stderr,
    )

    output_fq1_files, output_fq2_files = [mapping_fq1_files[0],], [mapping_fq2_files[0],]
    # Perform MASH distance filtering only when the number of samples is greater than 1.
    if len(mapping_fq1_files) > 1:
        print(
            "\n#############################MASH distance#####################################\n",
            file=sys.stderr,
        )
        # Preprocess: Calculate the MASH distance for each FASTQ against the one with the highest abundance, returning samples with dist < max_dist_threshold.
        preprocess_dir = outdir / "Preprocess"
        preprocess_dir.mkdir(exist_ok=True)
        
        output_fq1_files, output_fq2_files = filter_fq_by_mash(
            bin_id,
            mapping_fq1_files,
            mapping_fq2_files,
            preprocess_dir,
            max_abundance_index,
            max_dist_threshold,
            mash_path,
            mash_sketch_args,
            threads,
        )
        print(
            "\n#############################MASH distance#####################################\n",
            file=sys.stderr,
        )

    print(
        "\n#############################Merge output fastq#####################################\n",
        file=sys.stderr,
    )
    # Merge the qualifying FASTQs and output.
    prefix = str(outdir / f"{bin_id}_bwa_mash")
    output_fq1_file, output_fq2_file = cat_fastq_files(output_fq1_files, output_fq2_files, prefix)
    print(
        "\n#############################Merge output fastq#####################################\n",
        file=sys.stderr,
    )
    return output_fq1_file, output_fq2_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Generate the fastq file which can be used for bin reassembly.\n
        Example: 
            python preprocess_bin_assembly.py -i bin_id -a highest_abundance_fq1 -A highest_abundance_fq2 -1 sample1_fq1,..,sampleN_fq1 -2 sample1_fq2,..,sampleN_fq2 -b bin_fasta -r ref_fasta workdir
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--bin_id",
        type=str,
        help="Bin id.",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fqdir",
        type=str,
        help="The directory including all necessary clean fq.",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--highest_abundance_fq1",
        type=str,
        help="Path to the highest abundance fastq1.",
        required=True,
    )
    # parser.add_argument(
    #     "-A",
    #     "--highest_abundance_fq2",
    #     type=str,
    #     help="Path to the highest abundance fastq2.",
    #     required=True,
    # )
    parser.add_argument(
        "-1",
        "--samples_fastq1",
        type=str,
        help="Fastq1 files including the highest abundance sample, separated by comma.",
        required=True,
    )
    # parser.add_argument(
    #     "-2",
    #     "--samples_fastq2",
    #     type=str,
    #     help="Fastq2 files including the highest abundance sample, separated by comma.",
    #     required=True,
    # )
    parser.add_argument(
        "-b",
        "--bin_fasta",
        type=str,
        help="Path to the bin fasta file.",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--reference_fasta",
        type=str,
        nargs="+",
        help="Path to the reference fasta file, can be one or more.",
        required=True,
    )
    parser.add_argument(
        "workdir", 
        type=str, 
        help="Working directory."
    )
    parser.add_argument(
        "--use_single_sample",
        action="store_true",
        help="If set, only the highest abundance sample will be used.",
    )
    parser.add_argument(
        "--not_use_reference",
        action="store_true",
    )
    parser.add_argument(
        "--bwa_exe",
        type=str,
        #default=str(Path(__file__).parent.joinpath("bwa").resolve()),
        default=shutil.which("bwa"),
        help="Path to the bwa executable. Default: [%(default)s]",
    )
    parser.add_argument(
        "--bwa_mem_options",
        type=str,
        default="",
        help="Options for bwa mem mapping. Default: [%(default)s]",
    )
    parser.add_argument(
        "--samtools_exe",
        type=str,
        #default=str(Path(__file__).parent.joinpath("samtools").resolve()),
        default=shutil.which("samtools"),
        help="Path to the samtools executable. Default: [%(default)s]",
    )
    parser.add_argument(
        "--samtools_view_options",
        type=str,
        default="-bF 2316",
        help="Options for samtools view to extract mapping reads. This defaults representing filtering of unmapped, secondary and supplementary alignments.\n"
        "Default: [%(default)s]",
    )
    parser.add_argument(
        "--mash_exe",
        type=str,
        #default=str(Path(__file__).parent.joinpath("mash").resolve()),
        default=shutil.which("mash"),
        help="Path to the mash executable. Default: [%(default)s]",
    )
    parser.add_argument(
        "--mash_sketch_options",
        type=str,
        default="-s 1000000 -k 32 -m 2",
        help="Options for mash sketch. Default: [%(default)s]",
    )
    parser.add_argument(
        "--max_dist_threshold",
        type=float,
        default=0.2,
        help="Maximum distance threshold to filter unclosed samples. Default: [%(default)s]",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=32,
        help="Number of threads for each job. Default: [%(default)s]",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format=" %(asctime)s - %(levelname)s - %(message)s",
    )
    logging.debug("Start of preprocess bin reassembly.")
    start_time = time.time()

    bin_id = args.bin_id
    fq_dir = get_fq(args.fqdir)

    if args.highest_abundance_fq1 in fq_dir:
        fq_list = fq_dir[args.highest_abundance_fq1]
        max_abundance_fq1 = Path(fq_list[0]).resolve()
        max_abundance_fq2 = Path(fq_list[1]).resolve() if len(fq_list) > 1 else None

    else:
        raise ValueError(f"{args.highest_abundance_fq1} not exist in {args.fqdir}.")
    

    samples_fastq1 = []
    samples_fastq2 = []
    for sampleID in args.samples_fastq1.split(','):
        if sampleID in fq_dir:
            fq_list = fq_dir[args.highest_abundance_fq1]
            samples_fastq1.append(Path(fq_list[0]).resolve())
            samples_fastq2.append(Path(fq_list[1]).resolve() if len(fq_list) > 1 else None)

        else:
            raise ValueError(f"{sampleID} not exist in {args.fqdir}.")
    # max_abundance_fq1 = Path(args.highest_abundance_fq1).resolve()
    # max_abundance_fq2 = Path(args.highest_abundance_fq2).resolve()
    # samples_fastq1 = [Path(i).resolve() for i in args.samples_fastq1.split(',')]
    # samples_fastq2 = [Path(i).resolve() for i in args.samples_fastq2.split(',')]
    bin_fasta = Path(args.bin_fasta).resolve()
    ref_fastas = [Path(i).resolve() for i in args.reference_fasta]
    workdir = Path(args.workdir).resolve()
    bwa_path = Path(args.bwa_exe)
    bwa_mem_options = args.bwa_mem_options
    samtools_path = Path(args.samtools_exe)
    samtools_view_options = args.samtools_view_options
    mash_path = Path(args.mash_exe)
    mash_sketch_options = args.mash_sketch_options
    max_dist_threshold = args.max_dist_threshold
    threads = args.threads

    if len(samples_fastq1) != len(samples_fastq2):
        raise ValueError("The number of samples fastq1 and fastq2 must be the same.")

    workdir.mkdir(exist_ok=True)
    if len(samples_fastq1) != len(samples_fastq2):
        raise ValueError("samples_fastq1 and samples_fastq2 must be paired.")
    for i in range(len(samples_fastq1)):
        check_file_can_be_opened(samples_fastq1[i])
        check_file_can_be_opened(samples_fastq2[i])
    check_file_can_be_opened(bin_fasta)
    for ref_fasta in ref_fastas:
        check_file_can_be_opened(ref_fasta)
    check_exe_can_be_run(bwa_path)
    check_exe_can_be_run(samtools_path)
    check_exe_can_be_run(mash_path)

    logging.info(f"Bin ID: {bin_id}")
    logging.info(f"highest abundance sample fq1: {max_abundance_fq1}")
    logging.info(f"highest abundance sample fq2: {max_abundance_fq2}")
    logging.info(f"Bin fasta: {bin_fasta}")
    logging.info(f"reference fasta: {ref_fasta}")
    logging.info(f"work directory: {workdir}")
    logging.info(f"BWA executable: {bwa_path}")
    logging.info(f"BWA mem options: {bwa_mem_options}")
    logging.info(f"samtools executable: {samtools_path}")
    logging.info(f"samtools view options: {samtools_view_options}")
    logging.info(f"MASH executable: {mash_path}")
    logging.info(f"MASH sketch options: {mash_sketch_options}")
    logging.info(f"Maximum distance threshold: {max_dist_threshold}")

    # Remove the max abundance sample from the samples FASTQ, then input to the function.
    samples_fastq1.remove(max_abundance_fq1)
    samples_fastq2.remove(max_abundance_fq2)
    output_fq1_file, output_fq2_file = preprocess_sample_fastqs(
        bin_id,
        samples_fastq1,
        samples_fastq2,
        max_abundance_fq1,
        max_abundance_fq2,
        bin_fasta,
        ref_fastas,
        workdir,
        args.use_single_sample,
        args.not_use_reference,
        bwa_path,
        bwa_mem_options,
        samtools_path,
        samtools_view_options,
        mash_path,
        max_dist_threshold,
        mash_sketch_options,
        threads
    )

    end_time = time.time()
    logging.info(
        f"End of preprocess bin reassembly. Time elapsed: {(end_time - start_time) / 60:.5f} miniutes.\n\n"
    )