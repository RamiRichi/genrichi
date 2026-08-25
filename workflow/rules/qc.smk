# ── QC: adapter trimming + per-sample QC + aggregate MultiQC ─────────────────


rule fastp:
    input:
        r1=lambda wc: samples_df.loc[wc.sample, "fastq_r1"],
        r2=lambda wc: samples_df.loc[wc.sample, "fastq_r2"],
    output:
        r1="results/{sample}/qc/trimmed_R1.fastq.gz",
        r2="results/{sample}/qc/trimmed_R2.fastq.gz",
        html="results/{sample}/qc/fastp.html",
        json="results/{sample}/qc/fastp.json",
    threads: 4
    conda: "../envs/qc.yaml"
    log: "logs/{sample}/fastp.log"
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
        fastp=expand("results/{sample}/qc/fastp.json", sample=SAMPLES),
        flagstat=expand("results/{sample}/align/{sample}.flagstat", sample=SAMPLES),
        mosdepth=expand(
            "results/{sample}/align/{sample}.mosdepth.summary.txt", sample=SAMPLES
        ),
        markdup=expand(
            "results/{sample}/align/{sample}.markdup_metrics.txt", sample=SAMPLES
        ),
    output:
        html="results/multiqc/multiqc_report.html",
        data=directory("results/multiqc/multiqc_report_data"),
    conda: "../envs/qc.yaml"
    log: "logs/multiqc.log"
    params:
        search_dirs=lambda wc, input: " ".join(
            sorted({str(Path(f).parent.parent) for f in input.fastp})
        ),
    shell:
        """
        multiqc {params.search_dirs} \
            --outdir results/multiqc \
            --filename multiqc_report.html \
            --force \
            2>{log}
        """
