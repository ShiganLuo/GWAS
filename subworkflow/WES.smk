import os
module fastqc:
	snakefile: "../modules/fastqc/fastqc.smk"

cutadapt_config = {
        "indir": workdir,
        "outdir":  outdir
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
print(f"Cutadapt parameters: {cutadapt_config}")
use rule trimming_Paired from cutadapt as WES_trimming_Paired

bwa_mem2_confg = {
    "Procedure": {
        "bwaMem2": config.get("Procedure",{}).get("bwaMem2") or "bwa-mem2",
        "samtools": config.get("Procedure",{}).get("samtools") or "samtools",
    },
    "outdir": outdir,
    "fasta": fasta,
    "bwaMem2_index_prefix": config.get('bwaMem2_index_prefix') or None,
    "paired_samples": paired_samples,
    "single_samples": single_samples,

}
module bwa_mem2:
    snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
    config: bwa_mem2_confg
print(f"BWA-MEM2 parameters: {bwa_mem2_confg}")
use rule bwaMem2_index from bwa_mem2 as WES_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as WES_bwaMem2_alignment




