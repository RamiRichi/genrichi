# ── Annotation: VEP (COSMIC + ClinVar + gnomAD) → TSV table ──────────────────


rule vep_annotate:
    input:
        vcf="results/{sample}/calling/{sample}.pass.vcf.gz",
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
        cosmic_vcf=config["annotation"]["cosmic"]["vcf"],
        clinvar_vcf=config["annotation"]["clinvar"]["vcf"],
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf --compress_output bgzip \
            --stats_file {output.summary} \
            --cache --dir_cache {input.cache} \
            --assembly {params.genome} \
            --fasta {input.ref} \
            --fork {threads} \
            --af_gnomadg \
            --custom {params.cosmic_vcf},COSMIC,vcf,exact,0,CDS,AA,CNT \
            --custom {params.clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN \
            {params.extra} \
            2>{log}
        tabix -p vcf {output.vcf}
        """


rule vcf_to_table:
    input:
        vcf="results/{sample}/annotation/{sample}.vep.vcf.gz",
    output:
        tsv="results/{sample}/annotation/{sample}.variants.tsv",
    conda: "../envs/annotation.yaml"
    log: "logs/{sample}/vcf_to_table.log"
    script:
        "../scripts/vcf_to_table.py"
