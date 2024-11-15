<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-MetaflowX_logo_dark.png">
    <img alt="nf-core MetaflowX Logo" src="docs/images/nf-core-MetaflowX_logo_light.png" style="width: 75%; height: auto;">
  </picture>
</h1>


[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with slurm](https://img.shields.io/badge/run%20with-slurm-1AAEE8.svg?labelColor=000000)](https://www.schedmd.com)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14166585-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14166585)

# MetaflowX User Manual

**MetaflowX** is a Nextflow-based metagenomics analysis workflow that processes short-read sequences and/or assembled contigs to automatically generate taxonomic compositions, community functional profiles, non-redundant gene catalogs with functional annotations, and high-quality MAGs. This workflow is designed to simplify and accelerate metagenomic data analysis for researchers and bioinformaticians, offering a streamlined and scalable approach.

<p align="center">
    <img src="docs/images/workflow.png" alt="nf-core-metassembly workflow overview" width="100%"  height="auto">
</p>

## Contents

- [1. Pipeline Summary](#1-pipeline-summary)
- [2. Getting Startedy](#2-getting-started)
  - [2.1 Prerequisites](#21-prerequisites)
  - [2.2 Installation](#22-installation)
- [3. How to run](#3-how-to-run)
  - [3.1 Basic usage](#31-basic-usage)
  - [3.2 Demo runs](#32-demo-runs)
  - [3.3 Advance usage](#33-advance-usage)
- [4. Output](#4-output)
- [5. Support](#5-support)
- [6. Credits](#6-credits)
- [7. Citations](#7-citations)



## 1. Pipeline Summary

1. Quality control ( [`fastp`](https://github.com/OpenGene/fastp) [`Trimmomatic`](https://github.com/usadellab/Trimmomatic) [`Bowtie2`](https://github.com/BenLangmead/bowtie2))
2. Contig assembly ( [`SPAdes`](https://github.com/ablab/spades) [`MEGAHIT`](https://github.com/voutcn/megahit) )
3. Microbial taxonomy and metabolic function analysis ( [`MetaPhlAn`](https://github.com/biobakery/MetaPhlAn) [`HUMAnN`](https://github.com/biobakery/humann) [`Kraken2`](https://github.com/DerrickWood/kraken2) )
4. Gene catalog construction ( [`Prodigal`](https://github.com/hyattpd/Prodigal) [`CD-HIT`](https://github.com/weizhongli/cdhit) [`eggNOG-mapper`](https://github.com/eggnogdb/eggnog-mapper) [`antiSMASH`](https://github.com/antismash/antismash) [`BiG-MAP`](https://github.com/medema-group/BiG-MAP) )
5. Automated binning analysis ( [`MetaBAT2`](https://bitbucket.org/berkeleylab/metabat) [`CONCOCT`](https://github.com/BinPro/CONCOCT) [`SemiBin2`](https://github.com/BigDataBiology/SemiBin) [`MaxBin2`](https://sourceforge.net/projects/maxbin/) [`MetaBinner`](https://github.com/ziyewang/MetaBinner) [`COMEBin`](https://github.com/ziyewang/COMEBin) [`binny`](https://github.com/a-h-b/binny) [`DAS_Tool`](https://github.com/cmks/DAS_Tool) [`Checkm2`](https://github.com/chklovski/CheckM2) [`dRep`](https://github.com/MrOlm/drep) [`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk) [`CoverM`](https://github.com/wwood/CoverM) [`Deepurify`](https://github.com/ericcombiolab/Deepurify) [`COBRA`](https://github.com/linxingchen/cobra) )
6. Report generation ( [`Jinja`](https://github.com/pallets/jinja) [`MultiQC`](https://github.com/MultiQC/MultiQC) )


For the details implementation of each module in the pipeline, please refer to the [module description](docs/modules.md).


## 2. Getting Started

### 2.1 Prerequisites
To ensure smooth analysis with MetaflowX, we strongly recommend pre-building both the software components and the reference databases before starting your analysis.

- **Environment Setup**: Follow the [Environment Installation and Configuration](docs/dependencies.md) guide to create a Conda environment and manage dependencies.
- **Database Setup**: Refer to the [Database Installation and Configuration](docs/database.md) guide for downloading and configuring required databases.

> [!NOTE]
> If you are new to Nextflow and nf-core, check the [Nextflow installation guide](https://nf-co.re/docs/usage/installation). Ensure your setup passes the `-profile test` before processing real data.

  
### 2.2 Installation

1. **Clone the repository:**
  ```bash
   git clone https://github.com/01life/MetaflowX.git
  ```

## 3. How to run

### 3.1 Basic usage


1. Prepare a samplesheet `samplesheet.csv` with your input data that looks as follows:


```csv
id,raw_reads1,raw_reads2
S1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
S2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
```

2. Now, you can run the pipeline using:

```bash
nextflow run MetaflowX \
   -profile <docker/singularity/conda/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md). You can use the following command to see all the parameters of the pipeline.

   ```bash
   nextflow run MetaflowX --help
   ```


> [!NOTE]
> MetaflowX relies on plenty of tools and their databases. For detailed installation and configuration instructions, please refer to the [dependencies guide](docs/dependencies.md), [database guide](docs/database.md) and [version documentation](docs/version.md).


### 3.2 Demo runs
The `test` folder contains demo input file. These files consist of paired-end (PE) reads in FASTQ format. The following command is used to test mode1, which corresponds to the Quality Control (QC) function.

   ```bash
   nextflow run MetaflowX \
      -profile test \
      --outdir test_result 
   ```



### 3.3 Advance usage

MetaflowX supports both single-ended and paired data processing, with flexible execution patterns including `skip` and `mode` patterns. The platform features advanced capabilities for database modifications and parameter customization. For comprehensive tutorials and implementation guidelines, please refer to our [Usage Documentation](docs/usage.md).


## 4. Output

The results generated by MetaflowX include the following sections:

✤ **Quality control**
- 01.CleanData/

✤ **Contig assembly**
- 02.Contig

✤ **Microbial taxonomy and metabolic function analysis**
- 101.MetaPhlAn  
- 102.HUMAnN

✤ **Gene catalog construction**
- 03.Geneset  
- 04.GenesetProfile

✤ **Automated binning analysis**
- 05.BinSet  
- 06.BinsetProfile 

✤ **Report generation**
- 07.MultiQC
- MetaflowX_Report_*.html

✤ **Pipeline information**
- pipeline_info 

For more in-depth information on the Metaflow's results, refer to [Output Documentation](docs/output.md).

## 5. Support

Refer to the [MetaflowX tutorial](docs/usage.md) for an overview of analysis options and example runs.

See the [change log](CHANGELOG.md) for a detailed record of all updates, modifications, and improvements to MetaflowX.

For any questions, please visit the [MetaflowX GitHub page](https://github.com/01life/MetaflowX/issues).


## 6. Credits

MetaflowX was originally written by 👩‍💻Yang ying and 👩‍💻Liang lifeng.

We thank the following people for their extensive assistance in the development of this pipeline:

👨‍💻Long Shibin
👨Xie hailiang


<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## 7. Citations

If you use nf-core/mag for your analysis, please cite the preprint as follows:

> **MetaflowX: A Scalable and Resource-Efficient Workflow for Multi-Strategy Metagenomic Analysis**
>


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
