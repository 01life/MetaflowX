# Database Installation and Configuration

🚀 [MetaflowX User Manual](../README.md)

MetaflowX utilizes 15 reference databases crucial for taxonomic identification and functional annotation. These databases require an additional storage allocation of about 436.6 GB for pre-construction.

## Contents
 - [Making host genome index by Bowtie 2](#making-host-genome-index-by-bowtie-2)
 - [Downloading the MetaPhlAn database](#downloading-the-metaphlan-database-25gb)
 - [Downloading the HUMAnN database](#downloading-the-humann-database-30gb)
 - [Downloading the Kraken2 standard database](#downloading-the-kraken2-standard-database)
 - [Downloading the eggNOG-mapper database](#downloading-the-eggnog-mapper-database-20gb)
 - [Downloading the CheckM2 database](#downloading-the-checkm2-database-3gb)
 - [Download the GTDB-tk database](#download-the-gtdb-tk-database-110g)
 - [Downloading the BiG-MAP database](#downloading-the-big-map-database)
 - [Downloading the RGI database](#downloading-the-rgi-database)
 - [Downloading the VFDB database](#downloading-the-vfdb-database)
 - [Downloading the Deepurify database](#downloading-the-deepurify-database)

## Making host genome index by Bowtie 2

If you want to remove human reads from raw sequencing in the QC module, you will need to download and index the human genome.

First, lets download and merge the human genome hg38:

``` bash 
mkdir -p /path/to/hostDB
cd /path/to/hostDB
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

Using bowtie2-build, build the human genome index.

``` bash
mkdir -p hg38_bowtie2_index
bowtie2-build \
    --verbose \
    --threads 4 \
    hg38.fa.gz \
    hg38_bowtie2_index/<index_prefix>

# <index_prefix>: The prefix for the custom index filename; it is recommended to use `hg38` here.
```

Remember to set the `host_db` and `host_db_index` variables in the [config file](../nextflow.config)!

```bash
host_db="/path/to/hostDB/hg38_bowtie2_index"
host_db_index="<index_prefix>"
```


## Downloading the MetaPhlAn database (~25GB)
Here is an example using MetaPhlAn 4 with the database version mpa_vJun23_CHOCOPhlAnSGB_202307. Before installing the database, please ensure that MetaPhlAn is correctly installed and that you have activated the appropriate environment. Official tutorials can be found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4).

```bash
metaphlan --install --bowtie2db /path/to/MetaPhlAn
```

If you are unable to download the database automatically, you can download it manually.
```bash
echo mpa_vJun23_CHOCOPhlAnSGB_202307 > /mateflowx/db/MetaPhlAn/mpa_latest
nohup wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307_marker_info.txt.bz2 &
nohup wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307_species.txt.bz2 &
nohup wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.tar &
nohup wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJun23_CHOCOPhlAnSGB_202307_bt2.tar &

tar xvf mpa_vJun23_CHOCOPhlAnSGB_202307_bt2.tar
tar xvf mpa_vJun23_CHOCOPhlAnSGB_202307.tar
```

Remember to set the `mpa_db` and `mpa_index` variables in the [config file](../nextflow.config)! You can verify the mpa_index in the mpa_latest file.

```bash
mpa_db="/path/to/MetaPhlAn"
mpa_index="xxxxx"

# mpa_index: specify the id of the database version to use
```


## Downloading the HUMAnN database (~30GB)

Here is an example using HUMAnN 3 with the database version 201901b. Before installing the database, please ensure that MetaPhlAn is correctly installed and that you have activated the appropriate environment. Official tutorials can be found [here](https://github.com/biobakery/humann?tab=readme-ov-file#1-download-humann).

```bash
INSTALL_LOCATION=/path/to/humannDB/

#To download the ChocoPhlAn database (16.5GB):
humann_databases --download chocophlan full $INSTALL_LOCATION
#To download the full UniRef90 database (20.7GB, recommended):
humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION
#To download the EC-filtered UniRef90 database (0.9GB):
humann_databases --download uniref uniref90_ec_filtered_diamond $INSTALL_LOCATION
#To download the full UniRef50 database (6.9GB):
humann_databases --download uniref uniref50_diamond $INSTALL_LOCATION
#To download the EC-filtered UniRef50 database (0.3GB):
humann_databases --download uniref uniref50_ec_filtered_diamond $INSTALL_LOCATION
#To download the gene families into other functional categories
humann_databases --download utility_mapping full $INSTALL_LOCATION

parallel "gunzip -dc utility_mapping/{} | sed -e '1iID\tname' > utility_mapping/{.}" ::: map_ko_name.txt.gz map_pfam_name.txt.gz map_eggnog_name.txt.gz
parallel "gunzip -dc utility_mapping/{1} | sed -e '1iID\tname' > utility_mapping/{2}" ::: map_ec_name.txt.gz ::: map_level4ec_name.txt

humann_path=$(which humann)
map_metacyc=$(dirname ${humann_path})/../lib/python3.10/site-packages/humann/data/misc/map_metacyc-pwy_name.txt.gz
parallel "gunzip -dc {1} | sed -e '1iID\tname' > utility_mapping/{2}" ::: ${map_metacyc} ::: map_MetaCyc_name.txt
```

Remember to set the `humann_chocophlan_db`, `humann_protein_db` and `humann_map_db` variables in the [config file](../nextflow.config)!

```bash
humann_chocophlan_db="/path/to/humannDB/chocophlan"
humann_protein_db="/path/to/humannDB/uniprot"
humann_map_db="/path/to/humannDB/full_mapping"
```

## Downloading the Kraken2 standard database

You will need an estimated 120GB of RAM and ~128GB of space to download the Kraken2 database. Note that this is only needed if you intend on running the KRAKEN2.

``` bash
kraken2-build --standard --threads 24 --db /path/to/kraken2DB
```

Remember to set the `kraken2_db` variable in the [config file](../nextflow.config)!

```bash
kraken2_db="/path/to/kraken2DB"
```


## Downloading the eggNOG-mapper database (~20GB)
We create a Diamond index by default, encompassing both Bacteria and Archaea. Additionally, we download the eggNOG 5.0 database. To ensure proper functionality of eggNOG-mapper, please follow these instructions to place the index and database in two separate folders. Official tutorials can be found [here](https://github.com/eggnogdb/eggnog-mapper/wiki).

```bash
mkdir -p /path/to/eggnog/eggnog_5.0/
mkdir -p /path/to/eggnog/bact_arch/

cd /path/to/eggnog/bact_arch/
create_dbs.py -m diamond --dbname bact_arch --taxa Bacteria,Archaea --data_dir /path/to/eggnog/bact_arch/ -y
diamond makedb --in bact_arch.faa --db bact_arch.dmnd

cd /path/to/eggnog/eggnog_5.0/
download_eggnog_data.py --data_dir /path/to/eggnog/eggnog_5.0/ -y 
```

Remember to set the `eggnog_diamond_db` and `eggnog_mapper_db` variable in the [config file](../nextflow.config)!

```bash
eggnog_diamond_db="/path/to/eggnog/bact_arch/bact_arch.dmnd"
eggnog_mapper_db="/path/to/eggnog/eggnog_5.0"
```

## Downloading the CheckM2 database (~3GB)
Here is an example using CheckM2 with the version 1.0.1. Before installing the database, please ensure that CheckM2 is correctly installed and that you have activated the appropriate environment. Official tutorials can be found [here](https://github.com/chklovski/CheckM2).

``` bash
checkm2 database --download --path /path/to/CheckM2DB
```

Remember to set the `checkm2_db` variable in the [config file](../nextflow.config)!

```bash
checkm2_db="/path/to/CheckM2DB/CheckM2_database/uniref100.KO.1.dmnd"
```

## Download the GTDB-tk database (~110G)

Here is an example using GTDB-tk (version `2.1.1`) with the database version `r220`. After downloading and decompressing the large file `gtdbtk_data.tar.gz`, create a Mash reference sketch database (.msh). Official tutorials can be found [here](https://ecogenomics.github.io/GTDBTk/).

```bash
mkdir /path/to/gtdbtk_db
cd /path/to/gtdbtk_db

wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz

or

wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz

tar xvzf gtdbtk_data.tar.gz

perl -alne 'print "/path/to/gtdbtk_db/fastani/$F[1]$F[0]";' /path/to/gtdbtk_db/fastani/genome_paths.tsv > tmp.txt

mash sketch -l -p 16 tmp.txt -o gtdb_ref_sketch.msh -k 16 -s 5000

```

Remember to set the `mash_db` variable in the [config file](../nextflow.config)!

```bash
gtdbtk_db="/path/to/gtdbtk_db"
mash_db="/path/to/gtdbtk_db/mash"
```

MetaflowX supports using gtdb_to_ncbi_majority_vote.py to translate GTDB taxonomy to NCBI. However, you need to manually download bac120_metadata_r220.tsv.gz and ar53_metadata_r220.tsv.gz and place them in the metadata directory.
```bash
mkdir /path/to/gtdbtk_db/metadata/
cd /path/to/gtdbtk_db/metadata/
wget  https://data.gtdb.ecogenomic.org/releases/release220/220.0/ar53_metadata_r220.tsv.gz
wget  https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz

```
Remember to set the `gtdb_archaeal_metadata`  and `gtdb_bacterial_metadata` variable in the [config file](../nextflow.config)!

```bash
gtdb_archaeal_metadata="/path/to/gtdbtk_db/metadata/ar53_metadata_r220.tsv.gz"
gtdb_bacterial_metadata="/path/to/gtdbtk_db/metadata/bac120_metadata_r220.tsv.gz"
```

## Downloading the BiG-MAP database

To run BiG-SCAPE, you will also need to have the latest (processed) Pfam database Pfam-A.hmm.gz available from the Pfam FTP website (https://pfam.xfam.org/). Once the Pfam-A.hmm.gz file is downloaded, uncompress it and process it using the hmmpress command from the HMMER suit (http://hmmer.org/).


Remember to set the `bigspace_db` variable in the [config file](../nextflow.config)!

```bash
bigspace_db="/path/to/pfam_file"
```


## Downloading the RGI database

Download and build the latest AMR reference data from CARD to run RGI.

```bash
mkdir -p /path/to/CARD_db
cd /path/to/CARD_db
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json ./card.json 
```

Remember to set the `CARD_db` variable in the [config file](../nextflow.config)!

```bash
CARD_db="/path/to/CARD_db"
```


## Downloading the VFDB database

```bash
mkdir /path/to/VFDB_db
cd /path/to/VFDB_db
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
gzip -d VFDB_setB_pro.fas.gz
python simplify.seqID.py VFDB_setB_pro.fas VFDB_setB_pro.fas.S.fasta
diamond makedb --in VFDB_setB_pro.fas.S.fasta --db VFDB_setB_pro.fas.S
```

Remember to set the `VFDB_db` variable in the [config file](../nextflow.config)!

```bash
VFDB_db="/path/to/VFDB_db/VFDB_setB_pro.fas.S.fasta"
```

## Downloading the Deepurify database
```
cd /path/of/this/Deepurify-DB
wget https://drive.google.com/file/d/1FXpxoXFYHcX9QAFe7U6zfM8YjalxNLFk/view?usp=sharing
```
Remember to set the `deepurify_db` variable in the [config file](../nextflow.config)!

```bash
deepurify_db="/path/of/this/Deepurify-DB"
```