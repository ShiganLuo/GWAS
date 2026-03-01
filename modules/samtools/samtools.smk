outdir = config.get("outdir", "output")
onlykeep_properpair = config.get("onlykeep_properpair", False)

rule bam_flagstat:
    input:
        bam = outdir + "/bam/{genome}/{sample_id}.bam"
    output:
        flagstat = outdir + "/samtools/bam/flagstat/{genome}/{sample_id}.txt"
    log:
        outdir + "/log/samtools/{genome}/{sample_id}/flagstat.log"
    conda:
        "samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        {params.samtools} flagstat {input.bam} > {output.flagstat}
        """ 

rule bam_sort:
    input:
        bam = outdir + "/bam/{genome}/{sample_id}.bam"
    output:
        bam = temp(outdir + "/samtools/bam-sorted/{genome}/{sample_id}.bam"),
        bai = temp(outdir + "/samtools/bam-sorted/{genome}/{sample_id}.bam.bai")
    log:
        outdir + "/log/samtools/{genome}/{sample_id}/sort.log"
    threads: 4
    conda:
        "samtools.yaml"
    params:
        keep_proper_pair = "-f 2" if onlykeep_properpair else "",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"     
    shell:
        """
        {params.samtools} view -h {params.keep_proper_pair} -F 0x4  {input.bam} | {params.samtools} sort  -@ {threads}  > {output.bam} 2>{log}
        {params.samtools} index -@ {threads} {output.bam} 2>>{log}
        """

rule samtools_result:
    input:
        bam = outdir + "/samtools/bam-sorted/{genome}/{sample_id}.bam",
        bai = outdir + "/samtools/bam-sorted/{genome}/{sample_id}.bam.bai",
        flagstat = outdir + "/samtools/bam/flagstat/{genome}/{sample_id}.txt"

