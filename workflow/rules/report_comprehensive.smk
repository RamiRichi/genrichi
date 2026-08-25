# ── Comprehensive Report: SNV + CNV + MSI + TMB ────────────────────────────────


rule comprehensive_report:
    input:
        variants="results/{sample}/annotation/{sample}.somatic_variants.tsv",
        cnr="results/{sample}/cnv/{sample}.cnr",
        call_cns="results/{sample}/cnv/{sample}.call.cns",
        cnv_scatter="results/{sample}/cnv/{sample}-scatter.png",
        msi_score="results/{sample}/msi/{sample}.msi",
        tumor_fastp="results/{sample}/tumor/qc/fastp.json",
        normal_fastp="results/{sample}/normal/qc/fastp.json",
        tumor_flagstat="results/{sample}/tumor/align/{sample}.tumor.flagstat",
        normal_flagstat="results/{sample}/normal/align/{sample}.normal.flagstat",
        tumor_mosdepth="results/{sample}/tumor/align/{sample}.tumor.mosdepth.summary.txt",
        normal_mosdepth="results/{sample}/normal/align/{sample}.normal.mosdepth.summary.txt",
        tumor_markdup="results/{sample}/tumor/align/{sample}.tumor.markdup_metrics.txt",
        normal_markdup="results/{sample}/normal/align/{sample}.normal.markdup_metrics.txt",
    output:
        html="results/{sample}/report/{sample}_comprehensive_report.html",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/comprehensive_report.log"
    params:
        sample_id=lambda wc: wc.sample,
        patient_id=lambda wc: samples_df.loc[wc.sample, "patient_id"]
            if "patient_id" in samples_df.columns else wc.sample,
        sex=lambda wc: samples_df.loc[wc.sample, "sex"]
            if "sex" in samples_df.columns else "Unknown",
        tumor_type=lambda wc: samples_df.loc[wc.sample, "tumor_type"]
            if "tumor_type" in samples_df.columns else "Unknown",
        panel_name=config["panel"]["name"],
        company=config["report"]["company"],
        logo=config["report"]["logo"],
        show_synonymous=config["report"]["show_synonymous"],
        msi_threshold=config["msi"]["threshold"],
        tmb_coding_mb=config["tmb"]["coding_mb"],
        tmb_high_threshold=config["tmb"]["high_threshold"],
        cnv_amp_threshold=config["cnv"]["amp_threshold"],
        cnv_del_threshold=config["cnv"]["del_threshold"],
    script:
        "../scripts/generate_comprehensive_report.py"
