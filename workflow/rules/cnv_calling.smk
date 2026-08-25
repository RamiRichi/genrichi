# ── Copy Number Variation: mosdepth log2-ratio method ─────────────────────────
# Pure-Python approach — no cnvkit/htslib dependency needed.
# Uses mosdepth per-region depth files already produced by align_paired.smk.
#
# Algorithm (calculate_cnv.py):
#   1. Load per-region mean depth for tumor and normal (.regions.bed.gz)
#   2. Compute log2((tumor+1) / (normal+1)) per region
#   3. Rolling median smoothing (window=3)
#   4. Segment consecutive same-direction regions
#   5. Write CNR, CNS, call.cns (cnvkit-compatible TSV) + scatter PNG
#
# Runs in the pre-built report.yaml conda env (python=3.11, pandas, numpy, matplotlib).


rule calculate_cnv:
    """
    Panel-based CNV estimation from mosdepth per-region coverage.
    Requires only pandas + numpy + matplotlib — no cnvkit or htslib.
    """
    input:
        tumor_regions="results/{sample}/tumor/align/{sample}.tumor.regions.bed.gz",
        normal_regions="results/{sample}/normal/align/{sample}.normal.regions.bed.gz",
    output:
        cnr="results/{sample}/cnv/{sample}.cnr",
        cns="results/{sample}/cnv/{sample}.cns",
        call_cns="results/{sample}/cnv/{sample}.call.cns",
        scatter="results/{sample}/cnv/{sample}-scatter.png",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/cnv.log"
    params:
        amp_threshold=config["cnv"]["amp_threshold"],
        del_threshold=config["cnv"]["del_threshold"],
        min_probes=config["cnv"]["min_probes"],
    script:
        "../scripts/calculate_cnv.py"


rule cnv_export_bed:
    """Extract significant CNV calls (amp/del) to BED format for the report."""
    input:
        call_cns="results/{sample}/cnv/{sample}.call.cns",
    output:
        bed="results/{sample}/cnv/{sample}.cnv_calls.bed",
    conda: "../envs/report.yaml"
    log: "logs/{sample}/cnv_export.log"
    params:
        amp_thr=config["cnv"]["amp_threshold"],
        del_thr=config["cnv"]["del_threshold"],
        min_probes=config["cnv"]["min_probes"],
    run:
        import pandas as pd
        df = pd.read_csv(input.call_cns, sep="\t")
        # Keep only significant events with enough probes
        sig = df[
            ((df["type"] == "Amplification") | (df["type"] == "Deletion")) &
            (df["probes"] >= params.min_probes)
        ].copy()
        # Write BED: chrom, start, end, name (gene|type|log2), score (0), strand
        sig["name"] = (
            sig["gene"].astype(str) + "|" +
            sig["type"].astype(str) + "|" +
            sig["log2"].round(3).astype(str)
        )
        sig["score"] = 0
        sig["strand"] = "."
        bed = sig[["chromosome", "start", "end", "name", "score", "strand"]]
        bed.to_csv(output.bed, sep="\t", index=False, header=False)
        with open(log[0], "w") as fh:
            fh.write(f"Wrote {len(bed)} significant CNV calls to {output.bed}\n")
