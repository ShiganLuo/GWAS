# variants

为什么写这个流程，有现存nextflow的sarer非常完善，但是它有几个缺点:
- 不支持多分组样本比较，比如：有三个组别c1、e1、e2；要鉴定e1和e2的体细胞突变
- samplesheet的meta设计的不够简单

设计目标：

- 概念直观，上手成本低

- 规则组合清晰，认知负担小

- 依赖关系简洁，环境部署轻量


## 输入

meta_utils.py处理，多lane合并，单lane创建软链接

### 模式1

```python
python workflow/variants/utils/meta_utils.py --outdir output --fastq_dir data/Exome/Rawdata 
```
```
snakemake -s workflow/mutant-data-analyzer/main.smk --use-conda --cores 25 --config indir=data/Exome/Rawdata outdir=output metadata=data/Exome/samplesheet.csv --conda-prefix /data/pub/zhousha/env/mutation_0.1
```
- --conda-prefix固定pipeline的conda环境，以便下个项目复用


### 模式2

```python
python workflow/variants/utils/meta_utils.py --outdir output --meta data/Exome/samplesheet.csv
```
samplesheet.csv

| 列名        | 描述                              | 是否必须 |
| --------- | ------------------------------- | ---- |
| sample_id | 样本名                             | 是    |
| data_id   | 数据资源名            | 否    |
| fastq_1   | fastq file for read1            | 是    |
| fastq_2   | fastq file for read2            | 是    |
| design    | ctr_x or exp_x。相同x的ctr和exp将进行比较 | 是    |




## 项目结构

推荐的 Snakemake 友好目录（在现有基础上逐步迁移即可）：

```plain
.
├── main.smk            # 入口，include/use modules 或定义 subworkflow
├── input.yaml          # 用户输入，也可通过--config key=value……输入
├── main.py             # 可选：CLI 封装参数解析
├── modules/            # 所有可复用规则
├── subworkflows/       # 黑盒子流程（若仍需）
├── config/             # config.yaml/json, profiles/
├── envs/               # conda/mamba 环境
├── scripts/            # 规则脚本
├── utils/              # CLI 辅助，如 meta_utils
├── assets/             # 静态资源
├── docs/               # 文档
├── tests/              # 测试目录（可选）
```

说明：

- subworkflow通过module组合module下的smk中的规则来实现指定流程；通过config实现smk接口
- module每个smk所属目录包含三个文件：smk、yaml(conda环境)、schema.yaml(smk接口)
- main.smk引入metdata，定义生成目标文件函数来控制流程




## 参考

- [genome-seek](https://github.com/OpenOmics/genome-seek.git)