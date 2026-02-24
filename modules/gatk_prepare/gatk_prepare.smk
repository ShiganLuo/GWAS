from pathlib import Path

ASSAY = config.get("assay", "WES")
OUTDIR = Path(config.get("outdir", "output")) / ASSAY
REFERENCE = config.get("reference", {}).get("fasta") or config.get("genome")
INTERVAL = config.get("reference", {}).get("interval")
KNOWN_SITES = config.get("reference", {}).get("known_sites", [])
TMP_DIR = config.get("tmp_dir", "tmp")
GATK = config.get("Procedure", {}).get("gatk") or "gatk"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"


def get_bam_for_addReadsGroup(wildcards):
    return str(OUTDIR / "bam" / f"{wildcards.sample_id}.bam")


rule addReadsGroup:
    input:
        bam = get_bam_for_addReadsGroup
    output:
        bam = temp(lambda w: str(OUTDIR / "gatk_prepare" / "RG" / f"{w.sample_id}.bam")),
        bai = temp(lambda w: str(OUTDIR / "gatk_prepare" / "RG" / f"{w.sample_id}.bam.bai"))
    log:
        lambda w: str(OUTDIR / "log" / "gatk_prepare" / w.sample_id / "addReadsGroup.log")
    threads: 16
    conda:
        config["conda"]["gatk"]
    params:
        id = "{sample_id}",
        javaOptions = "--java-options -Xmx15G",
        RGLB = config.get("addReadsGroup", {}).get("RGLB") or "lib1",
        RGPL = config.get("addReadsGroup", {}).get("RGPL") or "illumina",
        RGPU = config.get("addReadsGroup", {}).get("RGPU") or "unit1",
        gatk = GATK,
        samtools = SAMTOOLS
    shell:
        """
        echo "sample_id: {wildcards.sample_id}" > {log}
        {params.gatk} AddOrReplaceReadGroups {params.javaOptions} \
            --INPUT {input.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB {params.RGLB} --RGPL {params.RGPL} --RGPU {params.RGPU} --RGSM {params.id} >> {log} 2>&1
        {params.samtools} index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule MarkDuplicates:
    input:
        bam = lambda w: str(OUTDIR / "gatk_prepare" / "RG" / f"{w.sample_id}.bam")
    output:
        bam = temp(lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}.bam")),
        bai = temp(lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}.bai")),
        metrics = temp(lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}_Markdup-metrics.txt"))
    log:
        lambda w: str(OUTDIR / "log" / "gatk_prepare" / w.sample_id / f"{w.sample_id}_MarkDuplicates.log")
    threads: 8
    conda:
        config["conda"]["gatk"]
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = GATK
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam}   \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --METRICS_FILE {output.metrics} > {log} 2>&1
        """


rule BaseRecalibrator:
    input:
        bam = lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}.bam"),
        bai = lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}.bai"),
        ref = REFERENCE,
        interval = INTERVAL
    output:
        table = lambda w: str(OUTDIR / "gatk_prepare" / "bqsr" / f"{w.sample_id}.recal_data.table")
    log:
        lambda w: str(OUTDIR / "log" / "gatk_prepare" / w.sample_id / "gatk_bqsr.log")
    threads: 8
    conda:
        config["conda"]["gatk"]
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = GATK,
        known_sites = " ".join([f"--known-sites {ks}" for ks in KNOWN_SITES])
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" BaseRecalibrator \
            -R {input.ref} --input {input.bam} \
            {params.known_sites} \
            {('-L ' + input.interval) if input.interval else ''} \
            -O {output.table}
        """


rule ApplyBQSR:
    input:
        bam = lambda w: str(OUTDIR / "gatk_prepare" / "bam-sorted-Markdup" / f"{w.sample_id}.bam"),
        table = lambda w: str(OUTDIR / "gatk_prepare" / "bqsr" / f"{w.sample_id}.recal_data.table"),
        ref = REFERENCE,
        interval = INTERVAL
    output:
        bam = lambda w: str(OUTDIR / "gatk_prepare" / "bqsr" / f"{w.sample_id}.sorted.markdup.BQSR.bam")
    log:
        lambda w: str(OUTDIR / "log" / "gatk_prepare" / w.sample_id / "ApplyBQSR.log")
    conda:
        config["conda"]["gatk"]
    threads:
        8
    params:
        javaOptions = f"--java-options -Xmx30G -Djava.io.tmpdir={TMP_DIR}",
        gatk = GATK,
        interval = INTERVAL
    shell:
        """
        {params.gatk} {params.javaOptions} ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            -bqsr {input.table} \
            {('-L ' + params.interval) if params.interval else ''} \
            -O {output.bam}
        """
