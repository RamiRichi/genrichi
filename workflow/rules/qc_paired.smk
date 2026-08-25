# ── QC: fastp for tumor AND normal reads ──────────────────────────────────────
# Uses {read_type} wildcard = "tumor" or "normal"
# Input FASTQs are looked up from samples_df columns: tumor_r1/r2, normal_r1/r2


rule fastp_paired:
    input:
        r1=lambda wc: samples_df.loc[wc.sample, f"{wc.read_type}_r1"],
        r2=lambda wc: samples_df.loc[wc.sample, f"{wc.read_type}_r2"],
    output:
        r1="results/{sample}/{read_type}/qc/trimmed_R1.fastq.gz",
        r2="results/{sample}/{read_type}/qc/trimmed_R2.fastq.gz",
        html="results/{sample}/{read_type}/qc/fastp.html",
        json="results/{sample}/{read_type}/qc/fastp.json",
    threads: 4
    conda: "../envs/qc.yaml"
    log: "logs/{sample}/{read_type}/fastp.log"
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            --length_required {config[trimming][fastp][min_length]} \
            --qualified_quality_phred {config[trimming][fastp][quality_threshold]} \
            --detect_adapter_for_pe \
            --correction \
            --overrepresentation_analysis \
            2>{log}
        """


rule multiqc:
    input:
        fastp=expand(
            "results/{sample}/{rtype}/qc/fastp.json",
            sample=SAMPLES, rtype=["tumor", "normal"],
        ),
        flagstat=expand(
            "results/{sample}/{rtype}/align/{sample}.{rtype}.flagstat",
            sample=SAMPLES, rtype=["tumor", "normal"],
        ),
        markdup=expand(
            "results/{sample}/{rtype}/align/{sample}.{rtype}.markdup_metrics.txt",
            sample=SAMPLES, rtype=["tumor", "normal"],
        ),
    output:
        html="results/multiqc/multiqc_report.html",
        data=directory("results/multiqc/multiqc_report_data"),
    conda: "../envs/qc.yaml"
    log: "logs/multiqc.log"
    params:
        search_dirs=lambda wc, input: " ".join(
            sorted({str(Path(f).parent.parent.parent) for f in input.fastp})
        ),
    shell:
        """
        multiqc {params.search_dirs} \
            --outdir results/multiqc \
            --filename multiqc_report.html \
            --exclude mosdepth \
            --force \
            2>{log}
        """
