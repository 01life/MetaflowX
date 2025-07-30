<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-MetaflowX_logo_dark.png">
    <img alt="nf-core MetaflowX Logo" src="images/nf-core-MetaflowX_logo_light.png" style="width: 75%; height: auto;">
  </picture>
</h1>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000\&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000\&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with slurm](https://img.shields.io/badge/run%20with-slurm-1AAEE8.svg?labelColor=000000)](https://www.schedmd.com)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14166585-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14166585)

# MetaflowX 用户手册

**MetaflowX** 是一个基于 Nextflow 构建的可扩展、模块化的宏基因组分析流程。支持短序列或组装后的 contig 输入，自动完成物种组成分析、功能注释、基因集构建及 MAG 恢复等关键分析。

<p align="center">
    <img src="images/workflow.png" alt="nf-core-metassembly workflow overview" width="100%" height="auto">
</p>

## 目录

* [1. 流程简介](#1-流程简介)
* [2. 快速开始](#2-快速开始)
  * [2.1 快速测试](#21-快速测试)
  * [2.2 环境要求](#22-环境要求)
  * [2.3 安装方法](#23-安装方法)
  * [2.4 环境与数据库检查](#24-环境与数据库安装检查)
* [3. 如何运行](#3-如何运行)
  * [3.1 基础用法](#31-基础用法)
  * [3.2 高级用法](#32-高级用法)
* [4. 输出内容](#4-输出内容)
* [5. 技术支持](#5-技术支持)
* [6. 贡献](#6-)
* [7. 引用信息](#7-引用信息)

---

## 1. 流程简介

MetaflowX 包含以下分析步骤：

1. **质控处理** ( [`fastp`](https://github.com/OpenGene/fastp) [`Trimmomatic`](https://github.com/usadellab/Trimmomatic) [`Bowtie2`](https://github.com/BenLangmead/bowtie2))
2. **Contig 组装** ( [`SPAdes`](https://github.com/ablab/spades) [`MEGAHIT`](https://github.com/voutcn/megahit) )
3. **基于比对的物种分类与代谢功能分析**( [`MetaPhlAn`](https://github.com/biobakery/MetaPhlAn) [`HUMAnN`](https://github.com/biobakery/humann) [`Kraken2`](https://github.com/DerrickWood/kraken2) )
4. **基因集构建** ( [`Prodigal`](https://github.com/hyattpd/Prodigal) [`CD-HIT`](https://github.com/weizhongli/cdhit) [`eggNOG-mapper`](https://github.com/eggnogdb/eggnog-mapper) [`antiSMASH`](https://github.com/antismash/antismash) [`BiG-MAP`](https://github.com/medema-group/BiG-MAP) )
5. **MAG 恢复与评估** ( [`MetaBAT2`](https://bitbucket.org/berkeleylab/metabat) [`CONCOCT`](https://github.com/BinPro/CONCOCT) [`SemiBin2`](https://github.com/BigDataBiology/SemiBin) [`MaxBin2`](https://sourceforge.net/projects/maxbin/) [`MetaBinner`](https://github.com/ziyewang/MetaBinner) [`COMEBin`](https://github.com/ziyewang/COMEBin) [`binny`](https://github.com/a-h-b/binny) [`MetaDecoder`](https://github.com/liu-congcong/MetaDecoder) [`Vamb`](https://github.com/RasmussenLab/vamb) [`DAS_Tool`](https://github.com/cmks/DAS_Tool) [`MAGScoT`](https://github.com/ikmb/MAGScoT) [`Checkm2`](https://github.com/chklovski/CheckM2) [`dRep`](https://github.com/MrOlm/drep) [`Galah`](https://github.com/wwood/galah?tab=readme-ov-file#galah) [`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk) [`CoverM`](https://github.com/wwood/CoverM) [`Deepurify`](https://github.com/ericcombiolab/Deepurify) [`COBRA`](https://github.com/linxingchen/cobra) )
6. **自动报告生成** ( [`Jinja`](https://github.com/pallets/jinja) [`MultiQC`](https://github.com/MultiQC/MultiQC) )

> 各模块详细介绍参见：[模块文档](docs/modules.md)

---

## 2. 快速开始

### 2.1 快速测试

如果系统已安装好 **Nextflow**，你可以通过以下命令快速体验流程功能：

1. 克隆代码仓库：

```bash
git clone https://github.com/01life/MetaflowX.git
```

2. 执行以下任一测试：

#### 1️⃣ 测试一：完整流程结构测试（使用 `stub` 模式，无需 Docker 或 Conda）

该测试执行完整流程结构但跳过实际计算，仅检查流程逻辑，**仅需 Nextflow**。

```bash
nextflow -bg run MetaflowX -stub -profile test_stub --outdir stub_remote > stub1.out
```

#### 2️⃣ 测试二：运行 `nf-core/fastp` 模块（需要 **Docker**）

运行内置的 `nf-core/fastp` 模块，需要 Docker 环境支持。

```bash
nextflow -bg run MetaflowX -profile test --outdir remote > remote.out
```

💡 两个测试都将在几分钟内完成，并在指定 `--outdir` 路径下生成日志和输出。

> 🚫 **注意**：上述测试仅用于验证流程功能，**不适用于真实生物学分析**。

---

### 2.2 环境要求

我们建议提前准备所有软件环境及数据库：

* 软件安装请参考：[环境依赖文档](docs/dependencies.md)
* 数据库下载与配置请参考：[数据库文档](docs/database.md)

---

### 2.3 安装方法

1. 克隆仓库：

```bash
git clone https://github.com/01life/MetaflowX.git
```

---

### 2.4 环境与数据库检查

> 💡 初学者请参考 [Nextflow 安装指南](https://nf-co.re/docs/usage/installation)，并优先运行 `-profile test` 测试确保配置正确。

流程内置了一组测试用的双端数据位于 `test/data/`，包括：

* `fullwork_test_sample1_reads1.fq.gz` 与 `fullwork_test_sample1_reads2.fq.gz`
* `fullwork_test_sample2_reads1.fq.gz` 与 `fullwork_test_sample2_reads2.fq.gz`

准备 `reads.csv` 输入文件（参考 [基础用法](#31-基础用法)），运行：

```bash
nextflow run MetaflowX \
   -profile <docker/singularity/conda/.../institute> \
   --input reads.csv \
   --outdir full_test
```

此步骤将：

* 检查工具是否可用
* 验证数据库路径是否配置正确
* 测试主要分析步骤能否正常执行

> ⚠️ **重点提醒**：
>
> 请务必正确配置以下文件：
>
> * `nextflow.config`：基础参数与数据库路径
> * `conf/modules.config`：各模块软件环境与参数
> * `conf/base.config`：计算资源（CPU、内存等）
>
> ✅ `-profile` 参数应根据本地环境选择，例如 `slurm`、`docker`、`local` 等。

---

## 3. 如何运行

### 3.1 基础用法

1. 准备输入样本表 `samplesheet.csv`，格式如下：

```csv
id,raw_reads1,raw_reads2
S1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
S2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
```

2. 执行流程：

```bash
nextflow run MetaflowX \
   -profile <docker/singularity/conda/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> ⚠️ **注意**
> 流程参数请通过命令行或 `-params-file` 提供；`-c` 提供的配置文件不能用于设置流程参数；
> 详见 [nf-core 配置文档](https://nf-co.re/usage/configuration#custom-configuration-files)

获取所有参数说明：

```bash
nextflow run MetaflowX --help
```

> 💡 MetaflowX 使用众多工具及数据库，具体依赖说明见：
> [环境依赖](docs/dependencies.md)、[数据库](docs/database.md)、[版本说明](docs/version.md)

---

### 3.2 高级用法

MetaflowX 支持以下功能：

* 支持单端或双端序列
* 通过 `--mode` 和 `--skip` 灵活选择模块执行
* 自定义数据库路径与参数配置

详见：[高级用法文档](docs/usage.md)

---

## 4. 输出内容

流程运行后会生成如下目录结构：

✤ **质控结果**

* 01.CleanData/

✤ **Contig 组装**

* 02.Contig

✤ **物种分类与功能分析**

* 101.MetaPhlAn
* 102.HUMAnN

✤ **基因集构建**

* 03.Geneset
* 04.GenesetProfile

✤ **MAG 自动分箱**

* 05.BinSet
* 06.BinsetProfile

✤ **报告生成**

* 07.MultiQC
* MetaflowX\_Report\_\*.html

✤ **流程信息**

* pipeline\_info

详见：[输出结果说明](docs/output.md)

---

## 5. 技术支持

* 查看 [MetaflowX 教程](docs/usage.md)
* 更新日志见 [CHANGELOG.md](CHANGELOG.md)
* 提交问题至 [GitHub Issues](https://github.com/01life/MetaflowX/issues)

---

## 6. 致谢

MetaflowX 由以下成员开发：

👩‍💻 杨颖
👩‍💻 梁丽凤

并获得以下成员的贡献与反馈：

👨‍💻 龙世斌
👨 谢海亮

---

## 7. 引用信息

如果你在研究中使用 MetaflowX，请引用：

> **MetaflowX: A Scalable and Resource-Efficient Workflow for Multi-Strategy Metagenomic Analysis**

更多引用信息请参见 [`CITATIONS.md`](CITATIONS.md)。
