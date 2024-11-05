# Environment Installation and Configuration

🚀 [MetaflowX User Manual](../README.md)

Due to the intricate nature of bioinformatics tools, creating a single operational environment is impractical. Instead, MetaflowX organizes these tools into several independent Conda environments, which requires approximately 6 GB of storage space.

For installing and managing the software in this pipeline, we recommend using [`conda`](https://conda.io/miniconda.html) or [`mamba`](https://github.com/mamba-org/mamba). Detailed installation instructions can be found in the following content.


## Contents
- [1. Basic tools](#1-basic-tools)
- [2. Specific tools](#2-specific-tools)
    - [2.1 CheckM2](#21-checkm2)
    - [2.2 antiSMASH](#22-antismash)
    - [2.3 RGI](#23-rgi)
    - [2.4 BiG-MAP](#24-big-map)
    - [2.5 MetaBinner](#25-metabinner)
    - [2.6 COMEBin](#26-comebin)
    - [2.7 binny](#27-binny)
    - [2.8 Deepurify](#28-deepurify)
- [3. Configuration](#3-configuration)

## 1. Basic tools

Most of the tools that the workflow depends on can be installed in the same conda environment. Refer to the command below to create a new conda environment from the [basic.yml](environment/basic.yml) configuration file.

``` bash
conda env create -f basic.yml -n MetaflowX
```


**Set the environment variable to the directory containing the GTDB-Tk reference data**

``` bash
conda env config vars set GTDBTK_DATA_PATH="/path/to/conda_env/gtdbtk_db" -n MetaflowX
```


## 2. Specific tools

Some tools rely on specific versions of packages, and they are installed in separate conda environments to avoid conflicts. You can install them using the yml file provided by us or refer to the official website of the tool for installation instructions. You can choose the way you like for installation.

### 2.1 CheckM2

Refer to the command below to create a new conda environment from the [checkm2.yml](environment/checkm2.yml) configuration file. The version of CheckM2 is 1.0.1 .

``` bash
conda env create -f checkm2.yml -n checkm2
source actiavte checkm2
pip install checkm2
```

For more detailes refer to the [official instructions](https://github.com/chklovski/CheckM2#installation).


### 2.2 antiSMASH

``` bash
conda create -n antismash -c  bioconda -c conda-forge  antismash
```

For more detailes refer to the [official instructions](https://docs.antismash.secondarymetabolites.org/install/).


### 2.3 RGI

``` bash
conda create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi
```

For more detailes refer to the [official instructions](https://github.com/arpcard/rgi#installation).


### 2.4 BiG-MAP
 
``` bash
cd /path/to/your/tools
git clone https://github.com/medema-group/BiG-MAP.git
conda env create -f BiG-MAP_process.yml BiG-MAP_process
chmod +x BiG-MAP/src/*.py
cp -f BiG-MAP/src/*.py /path/to/miniconda3/envs/BiG-MAP_process/bin
```

For more detailes refer to the [official instructions](https://github.com/medema-group/BiG-MAP#installation).


### 2.5 MetaBinner

```bash
conda create -n metabinner_env python=3.7.6
conda activate metabinner_env
conda install -c bioconda metabinner
```

For more detailes refer to the [official instructions](https://github.com/ziyewang/MetaBinner).


### 2.6 COMEBin

``` bash
conda create -n comebin_env
conda activate comebin_env
conda install -c conda-forge -c bioconda comebin
```

For more detailes refer to the [official instructions](https://github.com/ziyewang/COMEBin).


### 2.7 binny

``` bash
cd /path/to/software
git clone https://github.com/a-h-b/binny.git
cd binny
./binny -i config/config.init.yaml 
```

For more detailes refer to the [official instructions](https://github.com/a-h-b/binny).


If you need to use binny for binning , remember to set the `binny` and `binny_path` variables in the [config file](../nextflow.config)!

```bash
binny = true
binny_path = "/path/to/software/binny"
```


### 2.8 Deepurify
```
git clone https://github.com/ericcombiolab/Deepurify.git
conda env create -n deepurify -f deepurify-conda-env.yml

export CHECKM2DB="/path/to/database/uniref100.KO.1.dmnd"

conda activate deepurify
pip install Deepurify==2.3.7
```



## 3. Configuration

After installing the tools, remember to configure the absolute path of the conda environment. The basic tool environment is configured in the  [nextflow.config](../nextflow.config) file, while the specific tool environment is configured in the [modules.config](../conf/modules.config) file.

``` bash

#config the profiles
vi nextflow.config

#specifying the installation directory of conda environment globally
conda.enabled = true
process.conda = "/path/to/conda_env/MetaflowX"


#config the conda environment for specific tool
vi modules.config

withName: BIGMAP {
    conda = "/path/to/conda_env/bigmap"
}

withName: BINNY {
    conda = "/path/to/conda_env/binny"
}

withName: CHECKM2 {
    conda = "/path/to/conda_env/checkm2"
}

withName: RGI {
    conda = "/path/to/conda_env/rgi"
}

withName: METABINNER {
    conda = "/path/to/conda_env/metabinner"
}

withName: DASTOOL {
    conda = "/path/to/conda_env/dastool"
}

withName: ANTISMASH {
    conda = "/path/to/conda_env/antismash"
}

withName: GTDB {
    conda = "/path/to/conda_env/gtdbtk"
}

```

