# ── Alignment: bwa mem → sort → MarkDuplicates → BQSR (tumor AND normal) ──────
# Uses {read_type} wildcard = "tumor" or "normal"
# Read group SM tag is set to "{sample}_{read_type}" so Mutect2 can distinguish
# the two samples when both BAMs are provided.


rule bwa_mem_paired:
    input:
        r1="results/{sample}/{read_type}/qc/trimmed_R1.fastq.gz",
        r2="results/{sample}/{read_type}/qc/trimmed_R2.fastq.gz",
        ref=config["ref"]["genome"],
    output:
        bam=temp("results/{sample}/{read_type}/align/{sample}.{read_type}.raw.bam"),
    threads: 8
    conda: "../envs/align.yaml"
    log: "logs/{sample}/{read_type}/bwa_mem.log"
    params:
        rg=lambda wc: (
            f"@RG\\tID:{wc.sample}_{wc.read_type}"
            f"\\tSM:{wc.sample}_{wc.read_type}"
            f"\\tPL:{config['alignment']['platform']}"
            f"\\tLB:{wc.sample}_{wc.read_type}_lib1"
            f"\\tPU:{wc.sample}_{wc.read_type}"
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


rule mark_duplicates_paired:
    input:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.raw.bam",
    output:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.markdup.bam",
        bai="results/{sample}/{read_type}/align/{sample}.{read_type}.markdup.bam.bai",
        metrics="results/{sample}/{read_type}/align/{sample}.{read_type}.markdup_metrics.txt",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/{read_type}/mark_duplicates.log"
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


rule base_recalibrator_paired:
    input:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.markdup.bam",
        ref=config["ref"]["genome"],
        dbsnp=config["ref"]["dbsnp"],
        bed=config["panel"]["bed"],
    output:
        table="results/{sample}/{read_type}/align/{sample}.{read_type}.recal.table",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/{read_type}/base_recalibrator.log"
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


rule apply_bqsr_paired:
    input:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.markdup.bam",
        table="results/{sample}/{read_type}/align/{sample}.{read_type}.recal.table",
        ref=config["ref"]["genome"],
    output:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.final.bam",
        bai="results/{sample}/{read_type}/align/{sample}.{read_type}.final.bam.bai",
    conda: "../envs/align.yaml"
    log: "logs/{sample}/{read_type}/apply_bqsr.log"
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


rule flagstat_paired:
    input:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.final.bam",
    output:
        "results/{sample}/{read_type}/align/{sample}.{read_type}.flagstat",
    conda: "../envs/align.yaml"
    shell:
        "samtools flagstat {input.bam} > {output}"


rule mosdepth_paired:
    input:
        bam="results/{sample}/{read_type}/align/{sample}.{read_type}.final.bam",
        bai="results/{sample}/{read_type}/align/{sample}.{read_type}.final.bam.bai",
        bed=config["panel"]["bed"],
    output:
        summary="results/{sample}/{read_type}/align/{sample}.{read_type}.mosdepth.summary.txt",
        regions="results/{sample}/{read_type}/align/{sample}.{read_type}.regions.bed.gz",
    threads: 4
    conda: "../envs/align.yaml"
    log: "logs/{sample}/{read_type}/mosdepth.log"
    params:
        prefix="results/{sample}/{read_type}/align/{sample}.{read_type}",
    shell:
        """
        mosdepth \
            --threads {threads} \
            --by {input.bed} \
            {params.prefix} {input.bam} \
            2>{log}
        """
