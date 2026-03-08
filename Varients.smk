import os
SNAKEFILE_FULL_Varients = workflow.snakefile
SNAKEFILE_DIR_Varients = os.path.dirname(SNAKEFILE_FULL_Varients)
onlykeep_properpair = True
Varients = False
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_Varients ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
VarientsYaml = get_yaml_path("Varients")
configfile: VarientsYaml
logging.info(f"Include TEtranscripts config: {VarientsYaml}")

def check_tool_path_or_exit(tool_cmd, tool_name):
    """
    检查工具命令是否有效：
    1. 检查是否为空字符串。
    2. 如果是绝对路径，检查文件是否存在。
    """
    
    if not tool_cmd or not tool_cmd.strip():
        print(f"FATAL ERROR: Tool '{tool_name}' command is empty or not set in config.")
        sys.exit(1)
    
    # 如果路径包含目录分隔符，我们假设它是绝对路径或相对路径，需要检查存在性
    # 如果只是 'gatk' 或 'bwa-mem2' (即依赖系统PATH)，则跳过 os.path.exists 检查
    if os.sep in tool_cmd or os.altsep in tool_cmd:
        if not os.path.exists(tool_cmd):
            print(f"FATAL ERROR: Configured path for tool '{tool_name}' does not exist:")
            print(f"-> Path: {tool_cmd}")
            sys.exit(1)
        if os.path.isdir(tool_cmd):
            print(f"FATAL ERROR: Configured path for tool '{tool_name}' is a directory, not an executable file:")
            print(f"-> Path: {tool_cmd}")
            sys.exit(1)

    return True

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
    run:
        # ----------------------------------------------------
        # 步骤 1: 路径有效性检查 (使用 Python 检查)
        # ----------------------------------------------------
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        check_tool_path_or_exit(params.bwa_mem2_cmd, "BWA-MEM2")
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")

        # ----------------------------------------------------
        # 步骤 2: 执行命令 (使用 subprocess 或 shell)
        # ----------------------------------------------------

        # 注意：在 run 块中执行 shell 命令，需要使用 f-string 或 subprocess
        shell_commands = f"""
        # GATK CreateSequenceDictionary (.dict)
        echo "{params.gatk_cmd} CreateSequenceDictionary -R {input.fasta} -O {output.fa_dict}" >> {log}
        {params.gatk_cmd} CreateSequenceDictionary -R {input.fasta} -O {output.fa_dict} >> {log} 2>&1

        echo "{params.bwa_mem2_cmd} index {input.fasta}" >> {log}
        # BWA-MEM2 Index (.bwa-mem2.* 文件)
        {params.bwa_mem2_cmd} index {input.fasta} >> {log} 2>&1

        echo "{params.samtools_cmd} faidx {input.fasta}" >> {log}
        # samtools FASTA Index (.fai)
        {params.samtools_cmd} faidx {input.fasta} >> {log} 2>&1
        """
        
        # 使用 Snakemake 的内置 shell() 函数执行多行命令
        shell(shell_commands)

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
    run:
        check_tool_path_or_exit(params.bwa_mem2_cmd, "BWA-MEM2")
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")

        shell_command = f"""
        ({params.bwa_mem2_cmd} mem -T 0 -t {threads} {input.fasta} {params.reads_cmd} 2>&1) \
        | ({params.samtools_cmd} view -b - 2>&1) \
        > {output.bam} 
        >> {log} 2>&1
        """
        
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)

rule sort:
    input:
        bam = outdir + "/Varients/{genome}/bam/{sample_id}.bam"
    output:
        bam = temp(outdir + "/Varients/{genome}/bam-sorted/{sample_id}.bam"),
        bai = temp(outdir + "/Varients/{genome}/bam-sorted/{sample_id}.bam.bai")
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/sort.log"
    params:
        keep_proper_pair = config.get('params',{}).get('keep_proper_pair',False),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 4
    conda:
        config['conda']['run']
    run:
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")
        shell_command = f"""
        {params.samtools_cmd} view -h {params.keep_proper_pair} -F 0x4  {input.bam} | {params.samtools_cmd} sort -@ {threads}  > {output.bam} 2>{log}
        {params.samtools_cmd} index -@ {threads} {output.bam} 2>>{log}                
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)


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
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 10
    conda:
        config['conda']['run']        
    run:
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")
        shell_command = f"""
        {params.gatk_cmd} AddOrReplaceReadGroups {params.javaOptions} \
            --INPUT {input.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB cfDNA --RGPL BGI --RGPU DNBSEQ --RGSM {params.id} > {log.log} 2>&1
        {params.samtools_cmd} index -@ {threads} {output.bam} 2>> {log.log}
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)

# remove duplicate reads
rule dedup:
    input:
        bam = outdir + "/Varients/{genome}/RG/{sample_id}.bam"
    output:
        bam = temp(outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam"),
        bai = temp(outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam.bai"),
        metrics = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}_dedup-metrics.txt"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/MarkDuplicates.log"
    params:
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 16
    conda:
        config['conda']['run']
    run:
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")
        shell_command = f"""
        {params.gatk_cmd} MarkDuplicates \
            -I {input.bam} -O {output.bam} -M {output.metrics} \
            --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate > {log} 2>&1
        {params.samtools_cmd} index -@ {threads} {output.bam} 2>> {log}
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)


rule extract_exom_deduped:
    input:
        bam = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam",
        bed = lambda wildcards: config['genome'][wildcards.genome]['exon']['bed']
    output:
        exonBam = outdir + "/Varients/{genome}/exon/{sample_id}_exon.bam"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/extractExom.log"
    params:
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 5
    conda:
        config['conda']['run']
    run:
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")
        shell_command = f"""
            {params.samtools_cmd} view -b -L {input.bed} {input.bam} > {output.exonBam} 2>{log}
            {params.samtools_cmd} index -@ {threads} {output.exonBam} 2>> {log}
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)        

rule run_BQSR:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
        bam = outdir + "/Varients/{genome}/bam-sorted-deduped/{sample_id}.bam",
        knownVcf = lambda wildcards: config['genome'][wildcards.genome]['BQSR']['vcf']
    output:
        recal_table = outdir + "/Varients/{genome}/BQSR/{sample_id}_recal_data.table",
        BQSR_bam = outdir + "/Varients/{genome}/BQSR/{sample_id}_bqsr.bam"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/BQSR.log"
    params:
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk')
    threads: 5
    conda:
        config['conda']['run']
    run:
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        shell_command = f"""
            {params.gatk_cmd} BaseRecalibrator -I {input.bam} -R {input.fasta} --known-sites {input.knowVcf} -O {output.recal_table} > {log} 2>&1
            {params.gatk_cmd} {params.gatk_cmd} -I {input.bam} -R {input.fasta} --bqsr-recal-file {output.recal_table} -O {output.BQSR_bam}
            {params.samtools_cmd} index -@ {threads} {output.BQSR_bam} 2>> {log}
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)         

rule extract_exom_bqsr:
    input:
        bam = outdir + "/Varients/{genome}/BQSR/{sample_id}_bqsr.bam",
        bed = lambda wildcards: config['genome'][wildcards.genome]['exon']['bed']
    output:
        exonBam = outdir + "/Varients/{genome}/exon/{sample_id}_bqsr_exon.bam"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/extractExom.log"
    params:
        samtools_cmd = config.get('tools', {}).get('samtools', 'samtools')
    threads: 5
    conda:
        config['conda']['run']
    run:
        check_tool_path_or_exit(params.samtools_cmd, "Samtools")
        shell_command = f"""
            {params.samtools_cmd} view -b -L {input.bed} {input.bam} > {output.exonBam} 2>{log}
            {params.samtools_cmd} index -@ {threads} {output.exonBam} 2>> {log}
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)
 


def get_final_bam_path(wildcards):
    """为 HaplotypeCaller 动态生成输入 BAM 文件的完整路径"""
    perform_bqsr = config.get('params', {}).get('perform_bqsr', False)
    extract_exome = config.get('params', {}).get('extract_exome', False)
    if perform_bqsr and extract_exome:
        return f"{outdir}/Varients/{wildcards.genome}/exon/{wildcards.sample_id}_bqsr_exon.bam"
    elif perform_bqsr:
        return f"{outdir}/Varients/{wildcards.genome}/BQSR/{wildcards.sample_id}_bqsr.bam"
    elif extract_exome:
        return f"{outdir}/Varients/{wildcards.genome}/exon/{wildcards.sample_id}_exon.bam"
    else:
        return f"{outdir}/Varients/{wildcards.genome}/bam-sorted-deduped/{wildcards.sample_id}.bam"

rule HaplotypeCaller:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
        fai = lambda wildcards: input.fasta + ".fai",
        bam = get_final_bam_path,
        bai = lambda wildcards, input: input.bam + ".bai",
    output:
        vcf = outdir + "/Varients/{genome}/vcf/{sample_id}.vcf.gz"
    log:
        outdir + "/log/Varients/{genome}/{sample_id}/haplotypeCaller.log"
    conda:
        config['conda']['run']
    params:
        javaOptions = "--java-options -Xmx35G", # "--java-options -Xmx10G", # if omit this option, a default value will be used.
        gatk_cmd = config.get('tools', {}).get('gatk', 'gatk'),
        tmp_dir = config.get('tmp_dir', '/tmp'),
    threads: 25
    run:
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        shell_command = f"""
        {params.gatk_cmd} HaplotypeCaller {params.javaOptions} \
            -R {input.fasta} -I {input.bam} -O {output.vcf} \
            --tmp-dir {params.tmp_dir} > {log} 2>&1
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)        

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
    run:
        check_tool_path_or_exit(params.gatk_cmd, "GATK")
        shell_command = f"""
        {params.gatk_cmd} VariantFiltration {params.javaOptions}\
            -R {input.fasta} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS20 -filter "FS > 20.0" \
            --filter-name QD2 -filter "QD < 2.0" \
            --filter-name DP10 -filter "DP < 10.0" \
            --filter-name QUAL20 -filter "QUAL < 20.0" \
            > {log} 2>&1
        """
        shell(f"echo 'Running command: {shell_command}' >> {log}")
        shell(shell_command)    
