# ── ACMG 2015 Classification ──────────────────────────────────────────────────


rule acmg_classify:
    input:
        tsv="results/{sample}/annotation/{sample}.germline_variants.tsv",
    output:
        tsv="results/{sample}/annotation/{sample}.acmg_classified.tsv",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/acmg_classify.log"
    params:
        ba1_af=config["acmg"]["ba1_af"],
        bs1_af=config["acmg"]["bs1_af"],
        pm2_af=config["acmg"]["pm2_af"],
        lof_intolerant_genes=config["acmg"]["lof_intolerant_genes"],
    script:
        "../scripts/acmg_classifier.py"
