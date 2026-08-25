# ── Annotation: VEP (COSMIC optional + ClinVar + gnomAD) → TSV ────────────────
# COSMIC requires a registered download from cancer.sanger.ac.uk.
# If the file is absent the --custom flag is silently omitted so the
# pipeline continues without COSMIC IDs. Add the file later and rerun
# with -R vep_annotate_comprehensive to pick it up.

import os as _os


rule vep_annotate_comprehensive:
    input:
        vcf="results/{sample}/snv/{sample}.pass.vcf.gz",
        ref=config["ref"]["genome"],
        cache=config["annotation"]["vep"]["cache_dir"],
    output:
        vcf="results/{sample}/annotation/{sample}.vep.vcf.gz",
        tbi="results/{sample}/annotation/{sample}.vep.vcf.gz.tbi",
        summary="results/{sample}/annotation/{sample}.vep_summary.html",
    threads: 4
    conda: "../envs/annotation.yaml"
    log: "logs/{sample}/vep.log"
    params:
        genome=config["annotation"]["vep"]["genome_build"],
        extra=config["annotation"]["vep"]["extra"],
        clinvar_vcf=config["annotation"]["clinvar"]["vcf"],
        # Build COSMIC flag only if the file actually exists on disk
        cosmic_flag=lambda wc: (
            "--custom {vcf},COSMIC,vcf,exact,0,CDS,AA,CNT".format(
                vcf=config["annotation"]["cosmic"]["vcf"]
            )
            if config["annotation"]["cosmic"]["vcf"]
            and _os.path.isfile(config["annotation"]["cosmic"]["vcf"])
            else "# COSMIC VCF not found — skipping"
        ),
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --format vcf \
            --output_file {output.vcf} \
            --vcf --compress_output bgzip \
            --stats_file {output.summary} \
            --cache --dir_cache {input.cache} \
            --assembly {params.genome} \
            --fasta {input.ref} \
            --fork {threads} \
            --af_gnomadg \
            {params.cosmic_flag} \
            --custom {params.clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN \
            {params.extra} \
            2>{log}
        tabix -p vcf {output.vcf}
        """


rule vcf_to_somatic_table:
    input:
        vcf="results/{sample}/annotation/{sample}.vep.vcf.gz",
    output:
        tsv="results/{sample}/annotation/{sample}.somatic_variants.tsv",
    conda: "../envs/annotation.yaml"
    log: "logs/{sample}/vcf_to_somatic_table.log"
    script:
        "../scripts/vcf_to_somatic_table.py"
