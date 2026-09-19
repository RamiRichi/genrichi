# ── Runtime tool observation (Phase 5.3) ─────────────────────────────────────
# Records the tool versions actually in use, by running the tools themselves
# inside the same conda environments the analysis rules use. The report's
# provenance block reads these files; the annotation probe additionally fails
# the run if the VEP executable does not match the VEP cache release, before
# annotation starts. `priority` makes the scheduler run these tiny jobs first.

from runtime_probe import vep_cache_version_from_extra

_RUNTIME_PROBE = os.path.join(workflow.basedir, "scripts", "runtime_probe.py")


rule observe_tools_annotation:
    output:
        "results/runtime/annotation_tool_versions.json",
    conda: "../envs/annotation.yaml"
    log: "logs/runtime/observe_tools_annotation.log"
    priority: 100
    params:
        probe=_RUNTIME_PROBE,
        cache_dir=config["annotation"]["vep"]["cache_dir"],
        cache_version=lambda wc: vep_cache_version_from_extra(config["annotation"]["vep"].get("extra")) or "",
        assembly=config["annotation"]["vep"]["genome_build"],
    shell:
        """
        python "{params.probe}" --label annotation --tools vep,bcftools --out {output} \
            --vep-cache-dir "{params.cache_dir}" --vep-cache-version "{params.cache_version}" \
            --vep-assembly "{params.assembly}" 2>&1 | tee {log}
        """


rule observe_tools_alignment:
    output:
        "results/runtime/alignment_tool_versions.json",
    conda: "../envs/align.yaml"
    log: "logs/runtime/observe_tools_alignment.log"
    priority: 100
    params:
        probe=_RUNTIME_PROBE,
    shell:
        """
        python "{params.probe}" --label alignment --tools bwa,samtools,mosdepth --out {output} 2>&1 | tee {log}
        """


rule observe_tools_calling:
    output:
        "results/runtime/calling_tool_versions.json",
    conda: "../envs/calling.yaml"
    log: "logs/runtime/observe_tools_calling.log"
    priority: 100
    params:
        probe=_RUNTIME_PROBE,
    shell:
        """
        python "{params.probe}" --label calling --tools gatk,bcftools --out {output} 2>&1 | tee {log}
        """


rule observe_tools_qc:
    output:
        "results/runtime/qc_tool_versions.json",
    conda: "../envs/qc.yaml"
    log: "logs/runtime/observe_tools_qc.log"
    priority: 100
    params:
        probe=_RUNTIME_PROBE,
    shell:
        """
        python "{params.probe}" --label qc --tools fastp,multiqc --out {output} 2>&1 | tee {log}
        """
