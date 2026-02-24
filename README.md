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

项目结构参考nextflow

```plain
.
├── README.md 流程概要
├── assests 流程静态资源（meta示例……）
├── config 流程配置
├── docs 流程详细文档
├── main.py 启动脚本
├── rules 最小不可拆分规则（可以有多个规则组合，面板规则）
│   ├── local
│   └── smk
├── subworkflow 子流程
│   ├── local
│   └── smk
├── utils 流程工具类、函数
└── workflow 主流程
```

local代表自身开发，smk代表复用（暂无）




## 参考

- [genome-seek](https://github.com/OpenOmics/genome-seek.git)