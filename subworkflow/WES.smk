import os
module fastqc:
	snakefile: "../modules/fastqc/fastqc.smk"

cutadapt_params = {
        "indir": "output/raw_fastq",
        "outdir": "output"
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_params
print(f"Cutadapt parameters: {cutadapt_params}")
use rule trimming_Paired from cutadapt as WES_trimming_Paired



