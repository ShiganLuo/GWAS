import os
module fastqc:
	snakefile: "../modules/fastqc/fastqc.smk"

cutadapt_config = {
        "indir": workdir,
        "outdir":  outdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        }
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
print(f"Cutadapt parameters: {cutadapt_config}")
use rule trimming_Paired from cutadapt as WES_trimming_Paired

bwa_mem2_confg = {
    "Procedure": {
        "bwaMem2": config.get("Procedure",{}).get("bwaMem2"),
        "samtools": config.get("Procedure",{}).get("samtools")
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

samtools_config = {
    "Procedure": {
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "outdir": outdir
}

module samtools:
    snakefile: "../modules/samtools/samtools.smk"
    config: samtools_config
print(f"Samtools parameters: {samtools_config}")
use rule bam_sort from samtools as WES_bam_sort


gatk_prepare_config = {
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "addReadsGroup": {
        "RGLB": config.get("addReadsGroup", {}).get("RGLB"),
        "RGPL": config.get("addReadsGroup", {}).get("RGPL"),
        "RGPU": config.get("addReadsGroup", {}).get("RGPU")
    },
    "outdir": outdir,
    "fasta": fasta,
    "tmp_dir": config.get("tmp_dir")
}
module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
use rule addReadsGroup from gatk_prepare as WES_addReadsGroup
use rule MarkDuplicates from gatk_prepare as WES_MarkDuplicates

gatk_somatic_config = {
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
    "outdir": outdir,
    "fasta": fasta,
    "tmp_dir": config.get("tmp_dir"),
    "mutect2_parameters": config.get("mutect2_parameters")
}
module gatk_somatic:
    snakefile: "../modules/gatk/gatk_somatic/gatk_somatic.smk"
    config: gatk_prepare_config
use rule somaticMutect2 from gatk_somatic as WES_somaticMutect2
use rule gatk_index from gatk_prepare as WES_gatk_index


