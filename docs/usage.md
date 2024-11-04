# MetaflowX: Usage

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow -params-file option. Custom config files including those provided by the -c Nextflow option can be used to provide any configuration except for parameters, see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).


## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The pipeline supports both single-end and paired-end data. To ensure accurate parsing the input data, it is crucial that the column names in the CSV table correspond to the type of data being analyzed. Please refer to the following sections for more detailed information.


### Paired-end data

The pipeline analyzes paired-end data by default. You can specify a CSV samplesheet input file that contains the paths to your FASTQ files and additional metadata. The CSV file can support the following column names.

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id` | Custom sample ID. This entry will be identical for multiple sequencing libraries/runs from the same sample.  |
| `raw_reads1` | Full path to FastQ file for raw Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `raw_reads2` | Full path to FastQ file for raw Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `clean_reads1` | Full path to FastQ file for cleaned Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `clean_reads2` | Full path to FastQ file for cleaned Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `contig` | Full path to FastA file for assembled contigs. File has to have the extension ".fa". |

### Single-end data

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


## Running the pipeline 

MetaflowX has a crucial parameter `--mode` that controls the execution of the analysis. This parameter is prioritized above all others and defaults to 0, which means the entire workflow (all analysis modules) is executed by default. In this case, you can set some parameters to skip analysis modules that you don't need to execute, like `--skip_qc`, `--skip_marker` and `--skip_binning`. Note that these parameters only take effect when mode is set to 0. If you need to execute a specific analysis module, you can adjust the mode parameter accordingly. Please refer to the examples below for specific usage.


### Run the whole pipeline with default parameters `--mode 0`

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

#### Optional parameters
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

### Run quality control for microbiome analysis `--mode 1`

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

### Run flexible short read assembly `--mode 2`

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


### Run marker gene-based profiling for taxonomic composition and functional potential in microbiome analysis `--mode 3`

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

### Run unveiling community metabolism through gene prediction, annotation, and abundance estimation `--mode 4`

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

### Run streamlined strategies for MAG generation and abundance estimation in microbiome studies `--mode 5`

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

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [MetaflowX releases page](https://github.com/MetaflowX/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `slurm`
  - A generic configuration profile to be used with [Slurm](https://slurm.schedmd.com/overview.html/)
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### ⭐`-resume`

> **NB:** This is very useful useful to continue executions that was stopped by an error !


Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.
If you don't wanna use `-resume` and convince that the pipeline would not run error, you can just write 
```bash
process.scratch = true
```
into your config by using the paramter `-c` as follows. scratch will not output the middle files to your work dir.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
