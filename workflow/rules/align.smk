# ── Alignment: BWA-MEM2 → sort → MarkDuplicates → BQSR ───────────────────────


rule bwa_mem2:
    input:
        r1="results/{sample}/qc/trimmed_R1.fastq.gz",
        r2="results/{sample}/qc/trimmed_R2.fastq.gz",
        ref=config["ref"]["genome"],
    output:
        bam=temp("results/{sample}/align/{sample}.raw.bam"),
    threads: 8
    conda: "../envs/align.yaml"
    log: "logs/{sample}/bwa_mem2.log"
    params:
        rg=lambda wc: (
            f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}"
            f"\\tPL:{config['alignment']['platform']}"
            f"\\tLB:{wc.sample}_lib1\\tPU:{wc.sample}"
        ),
        extra=config["alignment"]["bwa_mem2"]["extra"],
    shell:
        """
        bwa mem \
            -t {threads} \
            -R '{params.rg}' \
            {params.extra} \
            {input.ref} {input.r1} {input.r2} \
        | samtools sort -@ {threads} -o {output.bam} - \
        2>{log}
        """


rule mark_duplicates:
    input:
        bam="results/{sample}/align/{sample}.raw.bam",
    output:
        bam="results/{sample}/align/{sample}.markdup.bam",
        bai="results/{sample}/align/{sample}.markdup.bam.bai",
        metrics="results/{sample}/align/{sample}.markdup_metrics.txt",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/mark_duplicates.log"
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            2>{log}
        mv {output.bam}.bai {output.bai} 2>/dev/null || true
        samtools index {output.bam} {output.bai}
        """


rule base_recalibrator:
    input:
        bam="results/{sample}/align/{sample}.markdup.bam",
        ref=config["ref"]["genome"],
        dbsnp=config["ref"]["dbsnp"],
        bed=config["panel"]["bed"],
    output:
        table="results/{sample}/align/{sample}.recal.table",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/base_recalibrator.log"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            -L {input.bed} \
            -O {output.table} \
            2>{log}
        """


rule apply_bqsr:
    input:
        bam="results/{sample}/align/{sample}.markdup.bam",
        table="results/{sample}/align/{sample}.recal.table",
        ref=config["ref"]["genome"],
    output:
        bam="results/{sample}/align/{sample}.final.bam",
        bai="results/{sample}/align/{sample}.final.bam.bai",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/apply_bqsr.log"
    shell:
        """
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} \
            2>{log}
        samtools index {output.bam} {output.bai}
        """


rule flagstat:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
    output:
        "results/{sample}/align/{sample}.flagstat",
    conda: "../envs/align.yaml"
    shell:
        "samtools flagstat {input.bam} > {output}"


rule mosdepth:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
        bai="results/{sample}/align/{sample}.final.bam.bai",
        bed=config["panel"]["bed"],
    output:
        summary="results/{sample}/align/{sample}.mosdepth.summary.txt",
        regions="results/{sample}/align/{sample}.regions.bed.gz",
    threads: 4
    conda: "../envs/align.yaml"
    log: "logs/{sample}/mosdepth.log"
    params:
        prefix="results/{sample}/align/{sample}",
    shell:
        """
        mosdepth \
            --threads {threads} \
            --by {input.bed} \
            --quantize 0:1:10:30:100: \
            {params.prefix} {input.bam} \
            2>{log}
        """
