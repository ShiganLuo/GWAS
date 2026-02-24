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
├── README.md
├── assets/
├── config/
│   └── profiles/
├── docs/
├── envs/                # 通用环境
├── modules/             # 可复用/可覆写的规则集合（替代当前 rules/local）
├── subworkflows/        # 需要整块调用的子流程（可选）
├── scripts/             # 规则脚本（若多，可放 modules/*/scripts）
├── utils/               # CLI 辅助，如 meta_utils
├── tests/               # 最小数据与 CI
├── Snakefile            # 入口：解析参数，选择 subworkflow 或 module 组合
└── main.py              # 可选封装 CLI，调用 snakemake -s Snakefile ...
```

说明：

- local 代表自研规则或子流程，smk 代表复用的外部 Snakemake 组件。
- 将入口 Snakefile 放在 workflow/ 方便 Snakedeploy、profile 与模块复用；保留 main.smk 便于渐进迁移。
- assets/、resources/、envs/、tests/ 可按需落地，先创建占位再逐步充实。




## 参考

- [genome-seek](https://github.com/OpenOmics/genome-seek.git)