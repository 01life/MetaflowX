#!/usr/bin/env python

import argparse
import logging
import os
from pathlib import Path
import subprocess
import shutil
import time
import threading
from typing import List, Dict, Set

import pandas as pd


class Task:
    def __init__(self, name: str, bash_file: Path, dependencies: List["Task"] = []):
        self.name = name
        self.bash_file = bash_file
        self.dependencies = dependencies
        self.job_id = None
        self.status = "Waiting"
        self.done_file = bash_file.with_suffix(".done")

    def is_finished(self):
        if self.status == "Done" or self.status == "Failed":
            return True
        return self.done_file.exists()

    def wait_for_finish(self):
        while not self.is_finished():
            time.sleep(3)

    def is_done(self):
        if self.status == "Done":
            return True
        return (
            self.done_file.exists()
            and self.done_file.read_text().strip() == "Mission accomplished!"
        )

    def is_failed(self):
        if self.status == "Failed":
            return True
        return (
            self.done_file.exists() and self.done_file.read_text().strip() == "Failed!"
        )

    def submit(self, queue_system_type: str):
        # Check if dependent tasks are completed.
        dependencies_has_failed = False
        for dependency in self.dependencies:
            while not dependency.is_finished():
                time.sleep(3)
            if dependency.is_failed():
                if dependency.status != "Failed":
                    logging.error(
                        f"Job {dependency.job_id} of {dependency.name} failed!"
                    )
                    dependency.status = "Failed"
                dependencies_has_failed = True
            elif dependency.is_done():
                if dependency.status != "Done":
                    dependency.status = "Done"
                    logging.info(f"Job {dependency.job_id} of {dependency.name} done!")
        # If any dependent tasks have failed, do not submit the task.
        if dependencies_has_failed:
            self.done_file.write_text("Failed!")
            self.status = "Failed"
            return
        # Submit the task.
        if queue_system_type == "slurm":
            submit_cmd = "sbatch"
        elif queue_system_type == "PBS" or queue_system_type == "pbs":
            submit_cmd = "qsub"
        else:
            raise NotImplementedError(f"Queue system {queue_system_type} is not supported!")
        result = subprocess.run(
            f"{submit_cmd} {self.bash_file}", shell=True, capture_output=True, text=True
        )
        output = result.stdout.strip()
        self.job_id = int(output.split()[-1])
        self.status = "Submit"
        logging.info(f"Submitted task {self.name} as job {self.job_id}")


class DrepParallel:
    def __init__(
        self,
        drep_exe: Path,
        genome_files: List[Path],
        workdir: Path,
        min_split_num: int,
        chunk_size: int,
        queue_system_type: str,
        projectname: str,
        queue: str,
        threads: int,
        max_mem: int,
        drep_options: str,
        restart_file: Path = None,
    ):
        self.drep_exe = drep_exe
        self.genome_files = genome_files
        self.workdir = workdir
        self.min_split_num = min_split_num
        self.chunk_size = chunk_size
        self.queue_system_type = queue_system_type
        self.projectname = projectname
        self.queue = queue
        self.threads = threads
        self.max_mem = max_mem
        self.drep_options = drep_options
        self.restart_file = restart_file

        self.tmpdir = self.workdir / "tmp"
        self.tmpdir.mkdir(exist_ok=True)
        self.shelldir = self.workdir / "shell"
        self.shelldir.mkdir(exist_ok=True)
        # Clear any leftover done files from the last run.
        for file in self.shelldir.glob("*.done"):
            file.unlink()
        self.tasks = []

        if self.restart_file:
            self._read_restart()
        else:
            self.restart_file = self.workdir / "restart.txt"
            self._create_commands()
            self._write_restart()

    def run(self):
        self._submit_tasks()
        total_tasks, done_tasks, fail_task = 0, 0, 0
        for roundtasks in self.tasks:
            for task in roundtasks:
                total_tasks += 1
                if task.status == "Done":
                    done_tasks += 1
                elif task.status == "Failed":
                    fail_task += 1
        logging.info(
            f"Total tasks: {total_tasks}, Done: {done_tasks}, Failed: {fail_task}"
        )
        if fail_task == 0:
            if len(self.genome_files) >= self.min_split_num:
                self._summary_split_result()
            else:
                self._summary_complete_result()
        else:
            raise RuntimeError("Some tasks failed!")

    def _read_restart(self):
        last_round = 0
        with open(self.restart_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                if not line:
                    continue
                round, bash_file, status, dependencies = line.split("\t")
                round = int(round)
                bash_file = Path(bash_file)
                if round != last_round:
                    self.tasks.append([])
                    last_round = round
                if dependencies != "None":
                    new_dependencies = []
                    for dependency in dependencies.split(","):
                        for task in self.tasks[-2]:
                            if str(task.bash_file) == dependency:
                                new_dependencies.append(task)
                                break
                    task = Task(bash_file.stem, bash_file, new_dependencies)
                else:
                    task = Task(bash_file.stem, bash_file)
                task.status = status
                self.tasks[-1].append(task)

    def _write_restart(self):
        with open(self.restart_file, "w") as o:
            o.write("#Round\tbash_file\tstatus\tdependencies\n")
            for i, roundtasks in enumerate(self.tasks):
                for task in roundtasks:
                    if len(task.dependencies) > 0:
                        dependencies = ",".join(
                            [
                                str(dependency.bash_file)
                                for dependency in task.dependencies
                            ]
                        )
                    else:
                        dependencies = "None"
                    o.write(f"{i+1}\t{task.bash_file}\t{task.status}\t{dependencies}\n")

    def _split_genome_files(self):
        self.genome_files = sorted(self.genome_files)
        genome_chunks = [
            self.genome_files[i : i + self.chunk_size]
            for i in range(0, len(self.genome_files), self.chunk_size)
        ]
        if len(genome_chunks[-1]) < self.chunk_size // 2:
            genome_chunks[-2].extend(genome_chunks[-1])
            genome_chunks.pop()
        return genome_chunks

    def _get_command_header(self, jobname: str):
        if self.queue_system_type == "slurm":
            return self._get_slurm_header(jobname)
        elif self.queue_system_type == "qsub":
            return self._get_qsub_header(jobname)
        elif self.queue_system_type == "sge":
            return self._get_sge_header(jobname)
        else:
            raise NotImplementedError(self.queue_system_type)

    def _get_slurm_header(self, jobname: str):
        outlog = self.shelldir / f"{jobname}.o"
        errlog = self.shelldir / f"{jobname}.e"
        return f"""#!/bin/bash
#SBATCH -J {jobname}
#SBATCH -A {self.projectname}
#SBATCH -p {self.queue}
#SBATCH -c {self.threads}
#SBATCH --mem={self.max_mem}G
#SBATCH -e {str(errlog)}
#SBATCH -o {str(outlog)}
"""
    
    def _get_qsub_header(self, jobname: str):
        outlog = self.shelldir / f"{jobname}.o"
        errlog = self.shelldir / f"{jobname}.e"
        return f"""#!/bin/bash
#PBS -N {jobname}
#PBS -A {self.projectname}
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.threads}
#PBS -l mem={self.max_mem}gb
#PBS -e {str(errlog)}
#PBS -o {str(outlog)}
"""
    
    def _get_sge_header(self, jobname: str):
        outlog = self.shelldir / f"{jobname}.o"
        errlog = self.shelldir / f"{jobname}.e"
        return f"""#!/bin/bash
#$ -N {jobname}
#$ -A {self.projectname}
#$ -q {self.queue}
#$ -pe smp {self.threads}
#$ -l mem={self.max_mem}gb
#$ -e {str(errlog)}
#$ -o {str(outlog)}
"""

    def _write_bash_file(self, bash_file: Path, header: str, *commands: List[str]):
        bash_file.touch()
        result_str = header + "\n"
        for command in commands:
            result_str += command + " && \\" + "\n"
        done_file = bash_file.with_suffix(".done")
        result_str += "\n" + "#" * 80 + "\n"
        result_str += (
            f'echo "Mission accomplished!" > {done_file} || \\\n'
            + f'echo "Failed!" > {done_file}\n'
        )
        with open(bash_file, "w") as o:
            o.write(result_str)

    def _create_commands(self):
        if len(self.genome_files) >= self.min_split_num:
            # Split the genome file into each group according to the specified size.
            drep1_lst = self._split_genome_files()
            round = 1
            self.tasks = [[], []]
            drep2_lst = []
            # Deduplicate within each group, submitting tasks to SLURM.
            for i, drep_genome in enumerate(drep1_lst):
                basename = f"dRep_r{round}_{i}"
                bash_file = self.shelldir / f"{basename}.sh"
                workdir = self.tmpdir / basename
                workdir.mkdir(exist_ok=True)
                header = self._get_command_header(basename)
                command = f'{self.drep_exe} dereplicate {workdir} -p {self.threads} --genomes {" ".join(map(str, drep_genome))} {self.drep_options}'
                self._write_bash_file(bash_file, header, command)
                self.tasks[0].append(Task(bash_file.stem, bash_file))
                drep2_lst.append(workdir.joinpath("dereplicated_genomes", f"*{drep_genome[0].suffix}"))
            # Perform pairwise comparisons for deduplication, submitting tasks to SLURM.
            round = 2
            for i in range(len(drep2_lst)):
                for j in range(i + 1, len(drep2_lst)):
                    basename = f"dRep_r{round}_{i}_{j}"
                    bash_file = self.shelldir / f"{basename}.sh"
                    workdir = self.tmpdir / basename
                    workdir.mkdir(exist_ok=True)
                    header = self._get_command_header(basename)
                    command = f"{self.drep_exe} dereplicate {workdir} -p {self.threads} --genomes {drep2_lst[i]} {drep2_lst[j]} {self.drep_options}"
                    self._write_bash_file(bash_file, header, command)
                    self.tasks[1].append(
                        Task(
                            bash_file.stem,
                            bash_file,
                            [self.tasks[0][i], self.tasks[0][j]],
                        )
                    )
        else:
            # Do not perform splitting when the number of input FASTA files is less than min_split_num; calculate directly.
            basename = "dRep"
            bash_file = self.shelldir / f"{basename}.sh"
            workdir = self.tmpdir / basename
            workdir.mkdir(exist_ok=True)
            header = self._get_command_header(basename)
            command = f'{self.drep_exe} dereplicate {workdir} -p {self.threads} --genomes {" ".join(map(str, self.genome_files))} {self.drep_options}'
            self._write_bash_file(bash_file, header, command)
            self.tasks.append([Task(bash_file.stem, bash_file)])

    def _remove_data_tables(self, data_tables_dir: Path):
        if data_tables_dir.exists():
            has_data_table = False
            for data_table in data_tables_dir.glob("*.csv"):
                data_table.unlink()
                has_data_table = True
            if has_data_table:
                logging.warning(f"Data tables of {data_tables_dir} exists, remove it!")

    def _submit_tasks(self):
        threads = []
        for roundtasks in self.tasks:
            for task in roundtasks:
                if task.status == "Done":
                    continue
                data_tables_dir = self.tmpdir / task.name / "data_tables"
                self._remove_data_tables(data_tables_dir)
                t = threading.Thread(target=task.submit(self.queue_system_type))
                t.start()
                threads.append(t)
            time.sleep(1)
        for thread in threads:
            thread.join()
            self._write_restart()

        # Wait for the completion of the last layer of tasks in the DAG.
        for wait_task in self.tasks[-1]:
            wait_task.wait_for_finish()
            if wait_task.is_done():
                wait_task.status = "Done"
                logging.info(f"Job {wait_task.job_id} of {wait_task.name} done!")
            elif wait_task.is_failed():
                wait_task.status = "Failed"
                logging.error(f"Job {wait_task.job_id} of {wait_task.name} failed!")
            self._write_restart()

    def _read_data_tables(self, data_tables_dir: Path) -> Dict[str, Set]:
        cluster_info = {}
        wdb = data_tables_dir / "Wdb.csv"
        cdb = data_tables_dir / "Cdb.csv"
        if not wdb.exists():
            raise FileNotFoundError(f"{wdb} not found!")
        if not cdb.exists():
            raise FileNotFoundError(f"{cdb} not found!")
        wdb_df = pd.read_csv(wdb, index_col=None, header=0)
        cdb_df = pd.read_csv(cdb, index_col=None, header=0)
        for _, row in wdb_df.iterrows():
            cluster_info.setdefault(row["genome"], set())
            cluster_id = row["cluster"]
            cluster_info[row["genome"]].update(
                cdb_df.loc[
                    (cdb_df["secondary_cluster"] == cluster_id)
                    & (cdb_df["genome"] != row["genome"]),
                    "genome",
                ].tolist()
            )
        return cluster_info

    def _summary_split_result(self):
        round1_tables = []
        for task in self.tasks[0]:
            data_tables_dir = self.tmpdir / task.name / "data_tables"
            cluster_info = self._read_data_tables(data_tables_dir)
            round1_tables.append(cluster_info)
        round2_tables = []
        # Perform pairwise comparisons on the deduplication results.
        for task in self.tasks[1]:
            data_tables_dir = self.tmpdir / task.name / "data_tables"
            cluster_info = self._read_data_tables(data_tables_dir)
            for i in task.name.split("_")[2:4]:
                i = int(i)
                for drep_genome1 in round1_tables[i]:
                    if drep_genome1 not in cluster_info:
                        for drep_genome2, genomes in cluster_info.items():
                            if drep_genome1 in genomes:
                                cluster_info[drep_genome2].update(
                                    round1_tables[i][drep_genome1]
                                )
                    else:
                        cluster_info[drep_genome1].update(
                            round1_tables[i][drep_genome1]
                        )
            round2_tables.append(cluster_info)
        # For dereplicated_genomes with the same ID, take the intersection, then take the union.
        # Generate cluster results based on the genomes results above.
        union_genomes = []
        result_cluster_info = {}
        intersect_ids_info = self._get_intersect_ids()
        for id in intersect_ids_info:
            tmp_intersect_ids = intersect_ids_info[id]
            intersect_genomes = set(round2_tables[tmp_intersect_ids[0]].keys())
            for tmp_intersect_id in tmp_intersect_ids[1:]:
                intersect_genomes &= set(round2_tables[tmp_intersect_id].keys())
            union_genomes.extend(intersect_genomes)
            for intersect_genome in intersect_genomes:
                result_cluster_info.setdefault(intersect_genome, set())
                for tmp_intersect_id in tmp_intersect_ids:
                    if intersect_genome in round2_tables[tmp_intersect_id]:
                        result_cluster_info[intersect_genome].update(
                            round2_tables[tmp_intersect_id][intersect_genome]
                        )
        #Write the cluster results to a file.
        result_cluster_file = self.workdir / "dRep_cluster.txt"
        with open(result_cluster_file, "w") as o:
            for drep_genome in result_cluster_info:
                result_str = (
                    drep_genome
                    + "\t"
                    + ";".join(result_cluster_info[drep_genome])
                    + "\n"
                )
                o.write(result_str)
        # Copy the deduplicated genomes to the result directory.
        drep_genome_dir = self.workdir / "dereplicated_genomes"
        if drep_genome_dir.exists():
            logging.warning(f"Result dir {drep_genome_dir} exists, remove it!")
            shutil.rmtree(drep_genome_dir)
        drep_genome_dir.mkdir()
        genome_path_info = {
            genome_path.name: genome_path for genome_path in self.genome_files
        }
        for drep_genome in sorted(union_genomes):
            shutil.copy(genome_path_info[drep_genome], drep_genome_dir)

    def _summary_complete_result(self):
        data_tables_dir = self.tmpdir / self.tasks[0][0].name / "data_tables"
        cluster_info = self._read_data_tables(data_tables_dir)
        # Write the cluster results to a file.
        result_cluster_file = self.workdir / "dRep_cluster.txt"
        with open(result_cluster_file, "w") as o:
            for drep_genome in cluster_info:
                result_str = (
                    drep_genome + "\t" + ";".join(cluster_info[drep_genome]) + "\n"
                )
                o.write(result_str)
        # Copy the deduplicated genomes to the result directory.
        drep_genome_dir = self.workdir / "dereplicated_genomes"
        if drep_genome_dir.exists():
            logging.warning(f"Result dir {drep_genome_dir} exists, remove it!")
            shutil.rmtree(drep_genome_dir)
        shutil.copytree(
            self.tmpdir / self.tasks[0][0].name / "dereplicated_genomes",
            drep_genome_dir,
        )

    def _get_intersect_ids(self) -> Dict[int, List[int]]:
        intersect_ids_info = {}
        for i, task in enumerate(self.tasks[1]):
            group_pair_ids = task.name.split("_")[2:4]
            intersect_ids_info.setdefault(group_pair_ids[0], [])
            intersect_ids_info.setdefault(group_pair_ids[1], [])
            intersect_ids_info[group_pair_ids[0]].append(i)
            intersect_ids_info[group_pair_ids[1]].append(i)
        return intersect_ids_info


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This script divide a big dRep clustering job into pieces and submit
        jobs to remote computers over a network to make it parallel.
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--drep_exe",
        type=str,
        default=str(Path(__file__).parent.joinpath("dRep").resolve()),
        help="Path to dRep executable. Default: [%(default)s]",
    )
    parser.add_argument(
        "-g",
        "--genomes",
        type=str,
        nargs="+",
        required=True,
        help="Genome files to be dereplicated.",
    )
    parser.add_argument("workdir", type=str, help="Working directory.")
    parser.add_argument(
        "-S", 
        "--min_split_num", 
        type=int, 
        default=800, 
        help="When the number of input fasta files is less than --min_split_num, splitting will not be performed."
        " Default: [%(default)s]",
    )
    parser.add_argument(
        "-s",
        "--split_size",
        type=int,
        default=200,
        help="Number of genomes in each piece."
        " Default: [%(default)s]",
    )
    parser.add_argument(
        "-T",
        "--queue_system_type",
        type=str,
        default="slurm",
        choices=["slurm", "PBS", "SGE"],
        help="Default: [%(default)s]",
    )
    parser.add_argument(
        "-p",
        "--projectname",
        type=str,
        default="NO22080020",
        help="Project name of the submit system. Default: [%(default)s]",
    )
    parser.add_argument(
        "-q",
        "--queue",
        type=str,
        default="q_32_64",
        help="Queue name of the submit system. Default: [%(default)s]",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=32,
        help="Number of threads for each job. Default: [%(default)s]",
    )
    parser.add_argument(
        "-m",
        "--max_mem",
        type=int,
        default=62,
        help="Max memory required for each job, Suffix with G. Default: [%(default)s]",
    )
    parser.add_argument(
        "--drep_options",
        type=str,
        default="-sa 0.95 -nc 0.3 --ignoreGenomeQuality",
        help="Other options for dRep dereplicate. Default: [%(default)s]",
    )
    parser.add_argument(
        "-r",
        "--restart_file",
        type=str,
        default=None,
        help="Restart file, used after a crash of run. Default: [%(default)s]",
    )
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(exist_ok=True)
    logging.basicConfig(
        filename=workdir / "dRep_parallel.log",
        level=logging.DEBUG,
        format=" %(asctime)s - %(levelname)s - %(message)s",
    )
    logging.debug("Start of dRep_para")
    start_time = time.time()

    drep_exe = Path(args.drep_exe).absolute()
    if not drep_exe.exists():
        raise FileNotFoundError(f"Cannot find dRep executable: {drep_exe}")
    if not os.access(drep_exe, os.X_OK):
        raise PermissionError(f"{drep_exe} is not executable.")

    genomes = []
    for genome in args.genomes:
        genome = Path(genome).absolute()
        if not genome.exists():
            raise FileNotFoundError(f"Cannot find genome file: {genome}")
        genomes.append(genome)

    restart_file = args.restart_file
    if restart_file is not None:
        restart_file = Path(restart_file).absolute()
        if not restart_file.exists():
            raise FileNotFoundError(f"Cannot find restart file: {restart_file}")

    drep = DrepParallel(
        drep_exe,
        genomes,
        workdir,
        args.min_split_num,
        args.split_size,
        args.queue_system_type,
        args.projectname,
        args.queue,
        args.threads,
        args.max_mem,
        args.drep_options,
        restart_file,
    )
    drep.run()

    end_time = time.time()
    logging.info(
        f"End of dRep_para. Time elapsed: {(end_time - start_time) / 60:.5f} miniutes."
    )
