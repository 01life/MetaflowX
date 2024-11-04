# nf-core/metassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.beta dev - [2023-03-31]

This accessible version of the pipeline💥

### `Added`

- skip_qc:跳过QC
- skip_assembly:跳过Assembly
- 镜像化: 新建metassembly队列，解决eggnogIO并发问题

### `Improved`

- Modulize: 支持从组装或基因集开始，支持跳过QC、组装、分箱。可以单独使用的模块：
   - 基因集+基因集注释 <<input:contigs>>
   - 其他：跑完对应任务后终止进程即可
   - 备注：暂不支持单独分箱，因为分箱和分箱注释必须要基因集结果输入。
- clean:删除了一些冗余注释
- Parameter Wiki: 新增参数文档可视化网页，但有效期仅为两周，需要持续更新

### `Fixed`
- Maxbin2: 将Maxbin2输出的.fasta后缀重命名为.fa，以便下一步checkm2正常运行
- Concoct: 输出列表列名重命名

### `Known Issues`
- bin丰度:是否需要计算bin丰度，有待增加

## v0.1 dev - [2023-03-39]

Initial release of nf-core/metassembly, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- 多管道合并: 由concat改为join，并行化

### `Fixed`
- wrong channel: BINNING 子流程中，contigs管道总会发生样本匹配错乱，即样本1会到样本2的工作目录下。且只有contig混乱。

### `Known Issues`
- EGGNOG: eggnog没有镜像时高IO延迟，非常慢。因此只使用001作为流程测试文件。正式版本要将.first()删除
- 脚本路径: scripts都在生产目录下，带绝对路径。理论上放到nf的bin/下可以使用相对路径，但目前仍未探索出使用方式。
- bin丰度:是否需要计算bin丰度，有待增加


