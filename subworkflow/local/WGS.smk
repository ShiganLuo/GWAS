from copy import deepcopy
from pathlib import Path

ASSAY = "WGS"
ASSAY_REF = config.get("reference", {}).get(ASSAY, {})

ASSAY_CONFIG = deepcopy(config)
ASSAY_CONFIG.update(
	{
		"assay": ASSAY,
		"outdir": config.get("outdir", "output"),
		"indir": FASTQ_PREPARED,
		"reference": {
			"fasta": ASSAY_REF.get("fasta"),
			"interval": ASSAY_REF.get("interval"),
			"known_sites": ASSAY_REF.get("known_sites", []),
			"gnomad": ASSAY_REF.get("gnomad"),
			"pon": ASSAY_REF.get("pon"),
			"bwa_index_prefix": ASSAY_REF.get("bwa_index_prefix"),
		},
		"tmp_dir": config.get("tmp_dir", "tmp"),
		"genome": ASSAY_REF.get("fasta"),
		"Procedure": config.get("procedure", {}),
	}
)


module cutadapt_wgs:
	snakefile: "../../rules/local/cutadapt/cutadapt.smk"
	config: ASSAY_CONFIG

module bwa_mem2_wgs:
	snakefile: "../../rules/local/bwa_mem2/bea_mem2.smk"
	config: ASSAY_CONFIG

module gatk_prepare_wgs:
	snakefile: "../../rules/local/gatk_prepare/gatk_prepare.smk"
	config: ASSAY_CONFIG

module gatk_germline_wgs:
	snakefile: "../../rules/local/gatk_germline.smk"
	config: ASSAY_CONFIG

module gatk_somatic_wgs:
	snakefile: "../../rules/local/gatk_somatic.smk"
	config: ASSAY_CONFIG


use rule trimming_Paired from cutadapt_wgs as wgs_trimming_Paired
use rule trimming_Single from cutadapt_wgs as wgs_trimming_Single
use rule bwaMem2_index from bwa_mem2_wgs as wgs_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2_wgs as wgs_bwaMem2_alignment
use rule addReadsGroup from gatk_prepare_wgs as wgs_addReadsGroup
use rule MarkDuplicates from gatk_prepare_wgs as wgs_MarkDuplicates
use rule BaseRecalibrator from gatk_prepare_wgs as wgs_BaseRecalibrator
use rule ApplyBQSR from gatk_prepare_wgs as wgs_ApplyBQSR
use rule HaplotypeCaller from gatk_germline_wgs as wgs_HaplotypeCaller
use rule filterHaplotypeCallerVcf from gatk_germline_wgs as wgs_filterHaplotypeCallerVcf
use rule somaticMutect2 from gatk_somatic_wgs as wgs_somaticMutect2
use rule filterMutect2Vcf from gatk_somatic_wgs as wgs_filterMutect2Vcf


ruleorder: wgs_filterMutect2Vcf > wgs_somaticMutect2 > wgs_filterHaplotypeCallerVcf > wgs_HaplotypeCaller > wgs_ApplyBQSR > wgs_BaseRecalibrator > wgs_MarkDuplicates > wgs_addReadsGroup > wgs_bwaMem2_alignment > wgs_bwaMem2_index > wgs_trimming_Paired > wgs_trimming_Single


rule wgs_prepare_fastq:
	input:
		meta = META
	output:
		prepared = directory(FASTQ_PREPARED)
	log:
		f"{OUTDIR}/WGS/log/prepare_fastq.log"
	params:
		fastq_dir = config.get("fastq", {}).get("source"),
		pipeline_outdir = OUTDIR
	shell:
		"""
		python {workflow.basedir}/../utils/meta_utils.py --meta {input.meta} --outdir {params.pipeline_outdir} --fastq_dir {params.fastq_dir} --log {log}
		"""


rule wgs_variants_core:
	input:
		germline = lambda w: f"{OUTDIR}/WGS/vcf-filtered/{w.sample}.vcf.gz",
		somatic = lambda w: f"{OUTDIR}/WGS/mutect2-vcf-filtered/{w.sample}.vcf.gz"
	output:
		marker = touch(lambda w: f"{OUTDIR}/WGS/variants/{w.sample}.done")
	log:
		lambda w: f"{OUTDIR}/WGS/log/{w.sample}/variants-core.log"
	params:
		assay = "WGS"
	shell:
		"""
		python {workflow.basedir}/../variants.py \
			--assay {params.assay} \
			--sample {wildcards.sample} \
			--germline_vcf {input.germline} \
			--somatic_vcf {input.somatic} \
			--log {log}
		touch {output.marker}
		"""
