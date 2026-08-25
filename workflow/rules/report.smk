# ── Report: aggregate QC + variants → per-sample HTML report ─────────────────


rule generate_report:
    input:
        variants="results/{sample}/annotation/{sample}.variants.tsv",
        fastp_json="results/{sample}/qc/fastp.json",
        flagstat="results/{sample}/align/{sample}.flagstat",
        mosdepth_summary="results/{sample}/align/{sample}.mosdepth.summary.txt",
        mosdepth_regions="results/{sample}/align/{sample}.mosdepth.regions.bed.gz",
        contamination="results/{sample}/calling/{sample}.contamination.table",
        markdup_metrics="results/{sample}/align/{sample}.markdup_metrics.txt",
    output:
        report="results/{sample}/report/{sample}_report.html",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/report.log"
    params:
        sample_id=lambda wc: wc.sample,
        panel_name=config["panel"]["name"],
        company=config["report"]["company"],
        min_vaf=config["report"]["min_vaf_report"],
        logo=config["report"]["logo"],
    script:
        "../scripts/generate_report.py"
