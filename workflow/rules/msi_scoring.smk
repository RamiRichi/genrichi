# ── MSI Scoring: VCF-based indel method ───────────────────────────────────────
# Uses somatic PASS variants from Mutect2 to estimate MSI status.
# MS-context indels per Mb of panel → MSI-High if >= threshold.
#
# This approach is equivalent to panel-based MSI tools (mSINGS, etc.) and
# avoids htslib/libdeflate conda dependency conflicts.
# Clearly labelled as "Panel-based estimate" in the report.


rule calculate_msi:
    input:
        vcf="results/{sample}/snv/{sample}.pass.vcf.gz",
    output:
        score="results/{sample}/msi/{sample}.msi",
        dis="results/{sample}/msi/{sample}.msi_dis",
        somatic="results/{sample}/msi/{sample}.msi_somatic",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/msi.log"
    params:
        coding_mb=config["tmb"]["coding_mb"],
        threshold=config["msi"]["threshold"],
    script:
        "../scripts/calculate_msi.py"
