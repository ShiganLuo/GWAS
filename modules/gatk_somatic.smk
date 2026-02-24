rule somaticMutect2:
    input:
        bam = outdir+"/bam-sorted-deduped/{sample_id}.bam",
        gn_fa_path=config['genome'],
        gnomad_path=outdir+"/snpRef/somatic-hg38_af-only-gnomad.hg38.vcf.gz",
        genome1k_path=outdir+"/snpRef/somatic-hg38_1000g_pon.hg38.vcf.gz"
    output:
        vcf=outdir+"/mutect2-vcf/{sample_id}.vcf.gz"
    log:
        outdir+"/log/{sample_id}/mutect2.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        # tmp="tmp",
        # java="" # "--java-options -Xmx10G", # if omit this option, a default value will be used.
    threads: 16
    shell:
        """
        gatk Mutect2 \
            -R {input.gn_fa_path} -I {input.bam} -O {output.vcf} \
            --germline-resource {input.gnomad_path} \
            --panel-of-normals {input.genome1k_path} \
            --native-pair-hmm-threads 10 \
            --tmp-dir {tmp_dir} > {log} 2>&1
        """

rule filterMutect2Vcf:
    input:
        vcf=outdir+"/mutect2-vcf/{sample_id}.vcf.gz",
        gn_fa_path=config['genome']
    output:
        vcf=outdir+"/mutect2-vcf-filtered/{sample_id}.vcf.gz",
    log:
        outdir+"/log/{sample_id}/mutect2-filtering.log"
    conda:
        config['conda']['cfDNA_base']
    threads: 16
    params:
        # java="--java-options -Xmx15G",
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.gn_fa_path} -V {input.vcf} -O {output.vcf} \
            -f-score-beta 1 \
            --tmp-dir {tmp_dir} > {log} 2>&1
        """ 