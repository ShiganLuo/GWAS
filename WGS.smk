import os
SNAKEFILE_FULL_WES = workflow.snakefile
SNAKEFILE_DIR_TEtranscripts = os.path.dirname(SNAKEFILE_FULL_WES)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_FULL_WES ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
WESYaml = get_yaml_path("WES")
configfile: WESYaml
logging.info(f"Include TEtranscripts config: {WESYaml}")

def get_output_path(input_fasta, extension):
    """生成带特定后缀的索引文件名"""
    return input_fasta + extension

rule Varients_index:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
    output:
        # GATK/Picard 索引 (.dict)
        fa_dict = lambda wildcards: get_output_path(input.fasta, ".dict"),
        # Samtools 索引 (.fai)
        fai = lambda wildcards: get_output_path(input.fasta, ".fai"),
        # BWA-MEM2 索引
        bwa_c2i = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.c2i"),
        bwa_idx = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.idx"),
        bwa_pac = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.pac"),
        bwa_amb = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.amb"),
        bwa_ann = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.ann")
    log:
        outdir + "/log/Varients/{genome}/Varients_index.log"
    params:
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        bwa_mem2_cmd = config.get('tools', {}).get('bwa_mem2', 'bwa-mem2'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    conda:
        config['conda']['run']
    shell:
        """
        # GATK CreateSequenceDictionary (.dict)
        echo "{params.gatk_cmd} CreateSequenceDictionary -R {input.fasta} -O {output.fa_dict}" > {log}
        {params.gatk_cmd} CreateSequenceDictionary -R {input.fasta} -O {output.fa_dict} >> {log} 2>&1

        echo "{params.bwa_mem2_cmd} index {input.fasta}" >> {log}
        # BWA-MEM2 Index (.bwa-mem2.* 文件)
        {params.bwa_mem2_cmd} index {input.fasta} >> {log} 2>&1

        echo "{params.samtools_cmd} faidx {input.fasta}" >> {log}
        # samtools FASTA Index (.fai)
        {params.samtools_cmd} faidx {input.fasta} >> {log} 2>&1
        """

def get_alignment_input(wildcards):
    """
    function: Dynamically determines the input file type: paired-end or single-end sequencing.
    Based on the paired_samples and single_samples lists.This function is called in the star_align rule.

    param: 
        wildcards: Snakemake wildcards object containing the sample_id.
        paired_samples = ['sample1', 'sample2', ...]
        single_samples = ['sample3', 'sample4', ...]
    These lists must be defined in the Snakefile or config file.

    return: A list of input file paths for the STAR alignment step. 
    """
    logging.info(f"[get_alignment_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        logging.info(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logging.info(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )

rule bwaMem2_alignment:
    input:
        reads = get_alignment_input,
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta'],
        bwa_c2i = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.c2i"),
        bwa_idx = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.idx"),
        bwa_pac = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.pac"),
        bwa_amb = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.amb"),
        bwa_ann = lambda wildcards: get_output_path(input.fasta, ".bwa-mem2.ann")
    output:
        bam = temp(outdir + "/Varients/{genome}/bam/{sample_id}.bam")
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/bwaMem2_alignment.log"
    params:
        reads_cmd = lambda wildcards, input: " ".join(input.reads),
        bwa_mem2_cmd = config.get('tools', {}).get('bwa_mem2', 'bwa-mem2'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 15
    conda:
        config['conda']['run']
    shell:
        """
        {params.bwa_mem2_cmd} mem -T 0 -t {threads} {input.fasta} {params.reads_cmd} 2> {log} \
                | {params.samtools_cmd} view -b - > {output.bam} 2>> {log}
        """

rule sort:
    input:
        bam = outdir + "/Varients/{genome}/bam/{sample_id}.bam"
    output:
        bam = temp(outdir + "/Varients/{genome}/bam-sorted/{sample_id}.bam"),
        bai = temp(outdir + "/Varients/{genome}/bam-sorted/{sample_id}.bam.bai")
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/sort.log"
    params:
        keep_proper_pair = "-f 2" if onlykeep_properpair else "",
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 4
    conda:
        config['conda']['run']
    shell:
        """
        {params.bwa_mem2_cmd} view -h {params.keep_proper_pair} -F 0x4  {input.bam} | {params.bwa_mem2_cmd} sort -@ {threads}  > {output.bam} 2>{log}
        {params.bwa_mem2_cmd} index -@ {threads} {output.bam} 2>>{log}
        """

rule addReadsGroup:
    input:
        bam = outdir + "/Varients/{genome}/bam-sorted/{sample_id}.bam"
    output:
        bam = temp(outdir + "/Varients/{genome}/RG/{sample_id}.bam"),
        bai = temp(outdir + "/Varients/{genome}/RG/{sample_id}.bam.bai")
    log:
        log = outdir + "/log/Varients/{genome}/{sample_id}/addReadsGroup.log"
    params:
        id = "{sample_id}",
        javaOptions = "--java-options -Xmx15G",
        tmp_dir = config["tmp_dir"],
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 10
    conda:
        config['conda']['run']        
    shell:
        """
        {params.gatk_cmd} AddOrReplaceReadGroups {params.javaOptions} \
            --TMP_DIR {params.tmp_dir} --INPUT {input.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB cfDNA --RGPL BGI --RGPU DNBSEQ --RGSM {params.id} > {log.log} 2>&1
        {params.samtools_cmd} index -@ {threads} {output.bam} 2>> {log.log}
        """
# remove duplicate reads
rule dedup:
    input:
        bam = outdir + "/Varients/{genome}/RG/{sample_id}.bam"
    output:
        bam = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam",
        bai = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam.bai",
        metrics = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}_dedup-metrics.txt"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/MarkDuplicates.log"
    params:
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 16
    conda:
        config['conda']['run']
    shell:
        """
        {params.gatk_cmd} MarkDuplicates \
            -I {input.bam} -O {output.bam} -M {output.metrics} \
            --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate > {log} 2>&1
        {params.samtools_cmd} index -@ {threads} {output.bam} 2>> {log}
        """
rule HaplotypeCaller:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
        bam = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam",
    output:
        vcf = outdir + "/Varients/{genome}/vcf/{sample_id}.vcf.gz"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/haplotypeCaller.log"
    conda:
        config['conda']['run']
    params:
        tmp="temp",
        javaOptions = "--java-options -Xmx35G", # "--java-options -Xmx10G", # if omit this option, a default value will be used.
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
    threads: 25
    shell:
        """
        {params.gatk_cmd} HaplotypeCaller {params.javaOptions} \
            -R {input.fasta} -I {input.bam} -O {output.vcf} \
            --tmp-dir {params.tmp} > {log} 2>&1
        """

rule filterHaplotypeCallerVcf:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta'],
        vcf = outdir + "/Varients/{genome}/vcf/{sample_id}.vcf.gz"
    output:
        vcf = outdir + "/Varients/{genome}/vcf-filtered/{sample_id}.vcf.gz"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/haplotypeCaller-filtering.log"
    conda:
        config['conda']['run']
    threads: 25
    params:
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        javaOptions = "--java-options -Xmx35G", # "--java-options -Xmx15G",
        tmp_dir = "temp"
    shell:
        """
        {params.gatk_cmd} VariantFiltration {params.javaOptions}\
            -R {input.fasta} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS20 -filter "FS > 20.0" \
            --filter-name QD2 -filter "QD < 2.0" \
            --filter-name DP10 -filter "DP < 10.0" \
            --filter-name QUAL20 -filter "QUAL < 20.0" \
            --tmp-dir {params.tmp_dir} > {log} 2>&1
        """ 
