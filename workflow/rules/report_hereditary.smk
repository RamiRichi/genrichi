# ── Hereditary HTML Report ────────────────────────────────────────────────────


rule hereditary_report:
    input:
        variants="results/{sample}/annotation/{sample}.acmg_classified.tsv",
        fastp_json="results/{sample}/qc/fastp.json",
        flagstat="results/{sample}/align/{sample}.flagstat",
        mosdepth_summary="results/{sample}/align/{sample}.mosdepth.summary.txt",
        mosdepth_regions="results/{sample}/align/{sample}.regions.bed.gz",
        markdup_metrics="results/{sample}/align/{sample}.markdup_metrics.txt",
    output:
        html="results/{sample}/report/{sample}_hereditary_report.html",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/hereditary_report.log"
    params:
        sample_id=lambda wc: wc.sample,
        panel_name=config["panel"]["name"],
        company=config["report"]["company"],
        logo=config["report"]["logo"],
        show_benign=config["report"]["show_benign"],
        patient_id=lambda wc: samples_df.loc[wc.sample, "patient_id"]
            if "patient_id" in samples_df.columns else wc.sample,
        sex=lambda wc: samples_df.loc[wc.sample, "sex"]
            if "sex" in samples_df.columns else "Unknown",
        indication=lambda wc: samples_df.loc[wc.sample, "indication"]
            if "indication" in samples_df.columns else "",
    script:
        "../scripts/generate_hereditary_report.py"
