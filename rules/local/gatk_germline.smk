rule HaplotypeCaller:
    input:
        bam = outdir + "/bam-sorted-deduped/{sample_id}.bam",
    output:
        vcf = outdir + "/vcf/{sample_id}.vcf.gz"
    log:
        outdir + "/log/{sample_id}/haplotypeCaller.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        gn_fa_path = config['genome'],
        tmp = "temp",
        java = "--java-options -Xmx35G" # "--java-options -Xmx10G", # if omit this option, a default value will be used.
    threads: 25
    shell:
        """
        gatk HaplotypeCaller {params.java} \
            -R {params.gn_fa_path} -I {input.bam} -O {output.vcf} \
            --tmp-dir {params.tmp} > {log} 2>&1
        """

rule filterHaplotypeCallerVcf:
    input:
        vcf = outdir + "/vcf/{sample_id}.vcf.gz"
    output:
        vcf = outdir + "/vcf-filtered/{sample_id}.vcf.gz"
    log:
        outdir + "/log/{sample_id}/haplotypeCaller-filtering.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 25
    params:
        gn_fa_path = config['genome'],
        java = "--java-options -Xmx35G", # "--java-options -Xmx15G",
        tmp_dir = "temp"
    shell:
        """
        gatk VariantFiltration {params.java}\
            -R {params.gn_fa_path} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS20 -filter "FS > 20.0" \
            --filter-name QD2 -filter "QD < 2.0" \
            --filter-name DP10 -filter "DP < 10.0" \
            --filter-name QUAL20 -filter "QUAL < 20.0" \
            --tmp-dir {params.tmp_dir} > {log} 2>&1
        """ 