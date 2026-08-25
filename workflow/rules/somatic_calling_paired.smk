# ── Somatic Variant Calling: Mutect2 (tumor-normal paired) ────────────────────


rule mutect2_paired:
    """Tumor-normal Mutect2 call. Normal SM tag = {sample}_normal."""
    input:
        tumor_bam="results/{sample}/tumor/align/{sample}.tumor.final.bam",
        tumor_bai="results/{sample}/tumor/align/{sample}.tumor.final.bam.bai",
        normal_bam="results/{sample}/normal/align/{sample}.normal.final.bam",
        normal_bai="results/{sample}/normal/align/{sample}.normal.final.bam.bai",
        ref=config["ref"]["genome"],
        bed=config["panel"]["bed"],
    output:
        vcf="results/{sample}/snv/{sample}.mutect2.vcf.gz",
        tbi="results/{sample}/snv/{sample}.mutect2.vcf.gz.tbi",
        stats="results/{sample}/snv/{sample}.mutect2.vcf.gz.stats",
        f1r2="results/{sample}/snv/{sample}.f1r2.tar.gz",
    threads: 4
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/mutect2.log"
    params:
        normal_name=lambda wc: f"{wc.sample}_normal",
        gnomad=lambda wc: (
            f"--germline-resource {config['ref']['gnomad']}"
            if config["ref"]["gnomad"] else ""
        ),
        pon=lambda wc: (
            f"--panel-of-normals {config['ref']['pon']}"
            if config["ref"]["pon"] else ""
        ),
        extra=config["calling"]["mutect2"]["extra"],
    shell:
        """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            --normal-sample {params.normal_name} \
            -O {output.vcf} \
            -L {input.bed} \
            --f1r2-tar-gz {output.f1r2} \
            {params.gnomad} \
            {params.pon} \
            {params.extra} \
            2>{log}
        """


rule learn_orientation_model_paired:
    input:
        f1r2="results/{sample}/snv/{sample}.f1r2.tar.gz",
    output:
        model="results/{sample}/snv/{sample}.read_orientation.tar.gz",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/learn_orientation.log"
    shell:
        """
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.model} \
            2>{log}
        """


rule get_pileup_tumor:
    input:
        bam="results/{sample}/tumor/align/{sample}.tumor.final.bam",
        gnomad=config["ref"]["gnomad"],
        bed=config["panel"]["bed"],
    output:
        table="results/{sample}/snv/{sample}.tumor.pileup.table",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/pileup_tumor.log"
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -L {input.bed} \
            -O {output.table} \
            2>{log}
        """


rule get_pileup_normal:
    input:
        bam="results/{sample}/normal/align/{sample}.normal.final.bam",
        gnomad=config["ref"]["gnomad"],
        bed=config["panel"]["bed"],
    output:
        table="results/{sample}/snv/{sample}.normal.pileup.table",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/pileup_normal.log"
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -L {input.bed} \
            -O {output.table} \
            2>{log}
        """


rule calculate_contamination_paired:
    input:
        tumor="results/{sample}/snv/{sample}.tumor.pileup.table",
        normal="results/{sample}/snv/{sample}.normal.pileup.table",
    output:
        table="results/{sample}/snv/{sample}.contamination.table",
        segments="results/{sample}/snv/{sample}.segments.table",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/contamination.log"
    shell:
        """
        gatk CalculateContamination \
            -I {input.tumor} \
            --matched-normal {input.normal} \
            --tumor-segmentation {output.segments} \
            -O {output.table} \
            2>{log}
        """


rule filter_mutect_calls_paired:
    input:
        vcf="results/{sample}/snv/{sample}.mutect2.vcf.gz",
        stats="results/{sample}/snv/{sample}.mutect2.vcf.gz.stats",
        orientation="results/{sample}/snv/{sample}.read_orientation.tar.gz",
        contamination="results/{sample}/snv/{sample}.contamination.table",
        segments="results/{sample}/snv/{sample}.segments.table",
        ref=config["ref"]["genome"],
    output:
        vcf="results/{sample}/snv/{sample}.filtered.vcf.gz",
        tbi="results/{sample}/snv/{sample}.filtered.vcf.gz.tbi",
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


rule pass_variants_paired:
    input:
        vcf="results/{sample}/snv/{sample}.filtered.vcf.gz",
    output:
        vcf="results/{sample}/snv/{sample}.pass.vcf.gz",
        tbi="results/{sample}/snv/{sample}.pass.vcf.gz.tbi",
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
