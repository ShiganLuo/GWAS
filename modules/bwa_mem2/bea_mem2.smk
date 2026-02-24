#Mapping with bwa
rule bwaMem2_index:
    input:
        fa = config['genome']
    output:
        index_prefix = outdir + "/{genome}/bwa_index/{genome}"
    log:
        outdir + "/log/{genome}/bwa_indexing.txt"
    threads: 15
    conda:
        config['conda']['align'] 
    params:
        bwa-mem2 = config.get("Procedure",{}).get("bwaMem2") or "bwa-mem2"
    shell:
        """
        {params.bwa-mem2} index -p {output.index_prefix} {input.fa} > {log} 2>&1
        """

rule bwaMem2_alignment:
    input:
        fastq1 = outdir + "/{genome}/cutadapt/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{genome}/cutadapt/{sample_id}_2.fq.gz",
        index_prefix = lambda wildcards: config['bwaMem2'][wildcards.genome]
    output:
        bam = outdir + "/{genome}/bam/{sample_id}.bam",
        sam = temp(outdir + "/{genome}/bam/{sample_id}.sam")
    log:
        outdir + "/log/{genome}/{sample_id}/bwa-alignment.txt"
    threads: 15
    conda:
        config['conda']['align'] 
    params:
        bwa-mem2 = config.get("Procedure",{}).get("bwaMem2") or "bwa-mem2",
        samtools = config.get("Procedure",{}).get("samtools") or "samtools"
    shell:
        """
        {params.bwa-mem2} mem \
        -T 0 \
        -t {threads} \
        {input.index_prefix} {input.fastq1} {input.fastq2} \
        -o {output.sam} \
        > {log} 2>&1

        {params.samtools} view -b {output.sam} > {output.bam} 2>{log}
        """
    