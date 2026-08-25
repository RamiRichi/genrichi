# ── Variant Calling: Mutect2 → contamination → orientation → FilterMutectCalls ─


rule mutect2:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
        ref=config["ref"]["genome"],
        bed=config["panel"]["bed"],
    output:
        vcf="results/{sample}/calling/{sample}.mutect2.vcf.gz",
        tbi="results/{sample}/calling/{sample}.mutect2.vcf.gz.tbi",
        stats="results/{sample}/calling/{sample}.mutect2.vcf.gz.stats",
        f1r2="results/{sample}/calling/{sample}.f1r2.tar.gz",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/mutect2.log"
    params:
        pon=lambda wc: (
            f"--panel-of-normals {config['ref']['pon']}"
            if config["ref"]["pon"]
            else ""
        ),
        gnomad=lambda wc: (
            f"--germline-resource {config['ref']['gnomad']}"
            if config["ref"]["gnomad"]
            else ""
        ),
        extra=config["calling"]["mutect2"]["extra"],
    shell:
        """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -L {input.bed} \
            {params.pon} \
            {params.gnomad} \
            --f1r2-tar-gz {output.f1r2} \
            {params.extra} \
            2>{log}
        """


rule learn_orientation_model:
    input:
        f1r2="results/{sample}/calling/{sample}.f1r2.tar.gz",
    output:
        model="results/{sample}/calling/{sample}.read_orientation.tar.gz",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/learn_orientation.log"
    shell:
        """
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.model} \
            2>{log}
        """


rule get_pileup_summaries:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
        gnomad=config["ref"]["gnomad"],
        bed=config["panel"]["bed"],
    output:
        table="results/{sample}/calling/{sample}.pileup_summaries.table",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/pileup_summaries.log"
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -L {input.bed} \
            -O {output.table} \
            2>{log}
        """


rule calculate_contamination:
    input:
        pileup="results/{sample}/calling/{sample}.pileup_summaries.table",
    output:
        table="results/{sample}/calling/{sample}.contamination.table",
        segments="results/{sample}/calling/{sample}.segments.table",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/contamination.log"
    shell:
        """
        gatk CalculateContamination \
            -I {input.pileup} \
            --tumor-segmentation {output.segments} \
            -O {output.table} \
            2>{log}
        """


rule filter_mutect_calls:
    input:
        vcf="results/{sample}/calling/{sample}.mutect2.vcf.gz",
        stats="results/{sample}/calling/{sample}.mutect2.vcf.gz.stats",
        orientation="results/{sample}/calling/{sample}.read_orientation.tar.gz",
        contamination="results/{sample}/calling/{sample}.contamination.table",
        segments="results/{sample}/calling/{sample}.segments.table",
        ref=config["ref"]["genome"],
    output:
        vcf="results/{sample}/calling/{sample}.filtered.vcf.gz",
        tbi="results/{sample}/calling/{sample}.filtered.vcf.gz.tbi",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/filter_mutect.log"
    params:
        min_af=config["calling"]["filter"]["min_af"],
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            --stats {input.stats} \
            --ob-priors {input.orientation} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --min-allele-fraction {params.min_af} \
            -O {output.vcf} \
            2>{log}
        """


rule pass_variants:
    input:
        vcf="results/{sample}/calling/{sample}.filtered.vcf.gz",
    output:
        vcf="results/{sample}/calling/{sample}.pass.vcf.gz",
        tbi="results/{sample}/calling/{sample}.pass.vcf.gz.tbi",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/pass_variants.log"
    params:
        min_depth=config["calling"]["filter"]["min_depth"],
        min_alt=config["calling"]["filter"]["min_alt_reads"],
    shell:
        """
        bcftools view \
            -f PASS \
            -e 'FORMAT/DP[0] < {params.min_depth} || FORMAT/AD[0:1] < {params.min_alt}' \
            {input.vcf} \
        | bcftools sort \
        | bgzip -c > {output.vcf} \
        2>{log}
        tabix -p vcf {output.vcf}
        """
