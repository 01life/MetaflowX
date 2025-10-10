# MetaflowX: Usage

🚀 [MetaflowX User Manual](../README.md)

## Contents
- [1. Pipeline parameters](#1-pipeline-parameters)
- [2. Samplesheet input](#2-samplesheet-input)
    - [2.1 Paired-end data](#21-paired-end-data)
    - [2.2 Single-end data](#22-single-end-data)
- [3. Running the pipeline](#3-running-the-pipeline)
    - [3.1 Standard workflow](#31-standard-workflow)
    - [3.2 'skip' pattern](#32-skip-pattern)
    - [3.3 'mode' pattern](#33-mode-pattern)
        - [3.3.1 Run quality control for microbiome analysis `--mode 1`](#331-run-quality-control-for-microbiome-analysis---mode-1)
        - [3.3.2 Run flexible short read assembly `--mode 2`](#332-run-flexible-short-read-assembly---mode-2)
        - [3.3.3 Run marker gene-based profiling for taxonomic composition and functional potential in microbiome analysis `--mode 3` ](#333-run-marker-gene-based-profiling-for-taxonomic-composition-and-functional-potential-in-microbiome-analysis---mode-3)
        - [3.3.4 Run unveiling community metabolism through gene prediction, annotation, and abundance estimation `--mode 4`](#334-run-unveiling-community-metabolism-through-gene-prediction-annotation-and-abundance-estimation---mode-4)
        - [3.3.5 Run streamlined strategies for MAG generation and abundance estimation in microbiome studies `--mode 5`](#335-run-streamlined-strategies-for-mag-generation-and-abundance-estimation-in-microbiome-studies---mode-5)
- [4. Core Nextflow arguments](#4-core-nextflow-arguments)
- [5. Custom configuration](#5-custom-configuration)
- [6. Azure Resource Requests](#6-azure-resource-requests)
- [7. Running in the background](#7-running-in-the-background)
- [8. Nextflow memory requirements](#8-nextflow-memory-requirements)



## 1. Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow -params-file option. Custom config files including those provided by the -c Nextflow option can be used to provide any configuration except for parameters, see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).


## 2. Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The pipeline supports both single-end and paired-end data. To ensure accurate parsing the input data, it is crucial that the column names in the CSV table correspond to the type of data being analyzed. Please refer to the following sections for more detailed information.


### 2.1 Paired-end data

The pipeline analyzes paired-end data by default. You can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata. The CSV file can support the following column names.

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id` | Custom sample ID. This entry will be identical for multiple sequencing libraries/runs from the same sample.  |
| `raw_reads1` | Full path to FastQ file for raw Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `raw_reads2` | Full path to FastQ file for raw Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `clean_reads1` | Full path to FastQ file for cleaned Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `clean_reads2` | Full path to FastQ file for cleaned Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `contig` | Full path to FastA file for assembled contigs. File has to have the extension ".fa". |

### 2.2 Single-end data

For single-end data analysis, you need to configure `--single_end` or `--single_end true` when running the pipeline. The input CSV file can support the following column names.

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id` | Custom sample ID. This entry will be identical for multiple sequencing libraries/runs from the same sample. |
| `raw_se` | Full path to FastQ file for raw single-end short read data. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `clean_se` | Full path to FastQ file for cleaned single-end short read data. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `contig` | Full path to FastA file for assembled contigs. File has to have the extension ".fa". |   

> #### Please note the following requirements:
> 
> - comma-seperated columns
> - Valid file extension: .csv
> - Sample IDs must be unique
> - FastQ files must be compressed (.fastq.gz, .fq.gz)
> - Column name is very strict
> - Within one samplesheet only cleaned reads and contig can be specified at the same time


## 3. Running the pipeline

MetaflowX has a crucial parameter `--mode` that controls the execution of the analysis. This parameter is prioritized above all others and defaults to 0, which means the entire workflow (all analysis modules) is executed by default. In this case, you can set some parameters to skip analysis modules that you don't need to execute, like `--skip_qc`, `--skip_marker` and `--skip_binning`. Note that these parameters only take effect when mode is set to 0. If you need to execute a specific analysis module, you can adjust the mode parameter accordingly. Please refer to the examples below for specific usage.


### 3.1 Standard workflow

Run the whole pipeline with default parameters `--mode 0`

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

#### Input for PE data

```console
id,raw_reads1,raw_reads2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

#### Input for SE data

```console
id,raw_se
sample3,sample3.fastq.gz
sample4,sample4.fastq.gz
```

### 3.2 'skip' pattern
- skip_qc 
- skip_marker
- skip_binning

To skip unnecessary analysis modules, use the above parameters according to your situation. These parameters can be used alone or in combination. Note that if you use the `--skip_qc` parameter, you need to provide the cleaned reads, and the column name of the CSV file need to be changed to: **clean_reads1** and **clean_reads2** or **clean_se**.


If you have cleaned paired-end data and need to perform assembly, gene set, and binning analysis. You can execute this pipeline with the following command:

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --skip_qc \
   --skip_marker
```

samplesheet.csv

```console
id,clean_reads1,clean_reads2
sample1,sample1_clean_R1.fastq.gz,sample1_clean_R2.fastq.gz
sample2,sample2_clean_R1.fastq.gz,sample2_clean_R2.fastq.gz
```
### 3.3 'mode' pattern

#### 3.3.1 Run quality control for microbiome analysis `--mode 1`

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 1 
```

#### Input for PE data

```console
id,raw_reads1,raw_reads2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

#### Input for SE data

```console
id,raw_se
sample3,sample3.fastq.gz
sample4,sample4.fastq.gz
```

#### Optional parameters
The pipeline provides two robust tools for the meticulous pre-processing of raw short sequencing reads, namely fastp (default setting) and Trimmomatic. If you need to use Trimmatic, refer to the following command to configure the `--qc_tool` parameter.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 1 \
   --qc_tool trimmomatic
```


<span id="mode2"></span>

#### 3.3.2 Run flexible short read assembly `--mode 2`

The assembly module supports both raw and clean reads. If raw reads are provided, the pipeline automatically performs quality control before assembly.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --mode 2
```

#### Input for PE data

```console
id,raw_reads1,raw_reads2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```
or

```console
id,clean_reads1,clean_reads2
sample1,sample1_clean_R1.fastq.gz,sample1_clean_R2.fastq.gz
sample2,sample2_clean_R1.fastq.gz,sample2_clean_R2.fastq.gz
```

#### Input for SE data

```console
id,raw_se
sample3,sample3.fastq.gz
sample4,sample4.fastq.gz
```
or 

```console
id,clean_se
sample3,sample3_clean.fastq.gz
sample4,sample4_clean.fastq.gz
```

#### Optional parameters
The pipeline supports metaSPAdes (default setting) and MEGAHIT to perform the assembly of a set of metagenomic reads. If you need to use MEGAHIT, refer to the following command to configure the `--assembly_tool` parameter.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 2 \
   --assembly_tool megahit
```


#### 3.3.3 Run marker gene-based profiling for taxonomic composition and functional potential in microbiome analysis `--mode 3`

This analysis module supports both raw and clean reads either. If raw reads are provided, the pipeline automatically performs quality control first. Please refer to [mode2](#mode2) for the format of the samplesheet.csv.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --mode 3
```

#### Optional parameters
The primary function of this module is to utilize MetaPhlAn and Kraken2 for the swift and accurate taxonomic profiling of microbial communities. Furthermore, Microbial functional profiles are determined using HUMAnN3. Please note that the HUMAnN3 analysis only works if the MetaPhlAn analysis is performed. The default seetings are as follows:

- metaphlan `true`
- humann `true`
- kraken `false`

If you only want to perform MetaPhlAn and Kraken2 and don't need the HUMAnN3 analysis, refer to the following command ：

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --mode 3
   --humann false
   --kraken 
```

<span id="mode4"></span>

#### 3.3.4 Run unveiling community metabolism through gene prediction, annotation, and abundance estimation `--mode 4`

This analysis module need to provide cleaned reads and assembled contig file, the pipeline will check whether the input are met based on the column name before analysis.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR>
   --mode 4
```

#### Input for PE data

```console
id,clean_reads1,clean_reads2,contig
sample1,sample1_clean_R1.fastq.gz,sample1_clean_R2.fastq.gz,sample1_contigs.fa
sample2,sample2_clean_R1.fastq.gz,sample2_clean_R2.fastq.gz,sample2_contigs.fa
```

#### Input for SE data

```console
id,clean_se,contig
sample3,sample3_clean.fastq.gz,sample3_contigs.fa
sample4,sample4_clean.fastq.gz,sample4_contigs.fa
```

#### Optional parameters

 Beyond general function annotation, MetaflowX supports special function annotation tools and databases, including the Virulence Factor Database (VFDB) and user-provided custom protein databases, annotated using DIAMOND blastp. Antibiotic resistance genes (ARG) are detected based on the Comprehensive Antibiotic Resistance Database (CARD) using Resistance Gene Identifier (RGI), and Biosynthetic Gene Clusters (BGCs) are identified using antiSMASH. A bowtie2-based pipeline performs a sequence search of CDS sequences against custom nucleic acid databases, such as DNA sequences of VFDB. You can configure the follwing parameters for special functional annotation.

- bigspace_db
- CARD_db
- VFDB_db
- nucleotide_db
- protein_db

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 4 \
   --bigspace_db /path/to/pfam_file \
   --CARD_db /path/to/CARD_db \
   --VFDB_db /path/to/VFDB_db/xxx.fasta \
   --nucleotide_db /path/to/custom_nt_db/xxx.fasta \
   --protein_db /path/to/custom_pro_db/xxx.fasta
```

#### 3.3.5 Run streamlined strategies for MAG generation and abundance estimation in microbiome studies `--mode 5`

This analysis module need to provide cleaned reads and assembled contig file either. Please refer to [mode4](#mode4) for the format of the samplesheet.csv.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 5
```

#### Optional parameters

The MetaflowX-Binning module includes several functions: contig binning, bins abundance estimation, bins taxonomy classification, bins function estimation, bins refine and bins reassembly. These functions configure parameters to manage the execution of their associated analysis. The default seetings are as follows:

- bin_abundance_calculation `true`
- bin_taxonomy `true`
- bin_function `true`
- bin_refine `false`
- bin_reassembly `false`

To enhance a set of bins by extracting reads specific to each bin and reassembling them, refer to the following command to configure the parameters.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 5 \
   --bin_refine \
   --bin_reassembly
```

It's essential to recognize that the quantity of bins is determined by the strategy of the binning algorithm. Consequently, MetaflowX offers seven different binning strategy binners: MetaBAT2, CONCOCT, SemiBin2, MaxBin2, MetaBinner, COMEBin, and Binny. The pipeline also sets parameters, and you can choose the most suitable binning software for analysis based on your actual data. The default seetings are as follows:

- metabat2 `true`
- concoct `true`
- semibin2 `true`
- maxbin2 `false`
- metabinner `false`
- binny `false`
- comebin `false`

If you need to use MetaBAT2, CONCOCT and MetaBinner for binning, refer to the following command to configure the parameters.

```bash
nextflow run MetaflowX \
   -profile slurm \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode 5 \
   --semibin2 false \
   --metabinner \
```
