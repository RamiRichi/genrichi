# ── Germline Calling: HaplotypeCaller (GVCF) → GenotypeGVCFs → Hard Filters ─


rule haplotypecaller:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
        bai="results/{sample}/align/{sample}.final.bam.bai",
        ref=config["ref"]["genome"],
        bed=config["panel"]["bed"],
        dbsnp=config["ref"]["dbsnp"],
    output:
        gvcf="results/{sample}/germline/{sample}.g.vcf.gz",
        tbi="results/{sample}/germline/{sample}.g.vcf.gz.tbi",
    threads: 4
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/haplotypecaller.log"
    params:
        extra=config["calling"]["haplotypecaller"]["extra"],
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -L {input.bed} \
            --dbsnp {input.dbsnp} \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            2>{log}
        """


rule genotype_gvcfs:
    input:
        gvcf="results/{sample}/germline/{sample}.g.vcf.gz",
        ref=config["ref"]["genome"],
        bed=config["panel"]["bed"],
    output:
        vcf="results/{sample}/germline/{sample}.genotyped.vcf.gz",
        tbi="results/{sample}/germline/{sample}.genotyped.vcf.gz.tbi",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/genotype_gvcfs.log"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            -L {input.bed} \
            2>{log}
        """


rule hard_filter_germline:
    """
    Apply GATK best-practice hard filters separately for SNPs and INDELs,
    then merge back and apply depth/GQ thresholds with bcftools.
    """
    input:
        vcf="results/{sample}/germline/{sample}.genotyped.vcf.gz",
        ref=config["ref"]["genome"],
    output:
        vcf="results/{sample}/germline/{sample}.filtered.vcf.gz",
        tbi="results/{sample}/germline/{sample}.filtered.vcf.gz.tbi",
    conda: "../envs/calling.yaml"
    log: "logs/{sample}/hard_filter_germline.log"
    params:
        min_dp=config["calling"]["filter"]["min_depth"],
        min_gq=config["calling"]["filter"]["min_gq"],
        tmpdir=lambda wc: f"results/{wc.sample}/germline/tmp",
    shell:
        """
        mkdir -p {params.tmpdir}

        # --- SNP hard filter ---
        gatk SelectVariants \
            -R {input.ref} -V {input.vcf} \
            --select-type-to-include SNP \
            -O {params.tmpdir}/snps.vcf.gz 2>{log}

        gatk VariantFiltration \
            -R {input.ref} \
            -V {params.tmpdir}/snps.vcf.gz \
            --filter-expression "QD < 2.0"             --filter-name "QD2"             \
            --filter-expression "FS > 60.0"            --filter-name "FS60"            \
            --filter-expression "MQ < 40.0"            --filter-name "MQ40"            \
            --filter-expression "MQRankSum < -12.5"    --filter-name "MQRankSum-12.5"  \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {params.tmpdir}/snps_filt.vcf.gz 2>>{log}

        # --- INDEL hard filter ---
        gatk SelectVariants \
            -R {input.ref} -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {params.tmpdir}/indels.vcf.gz 2>>{log}

        gatk VariantFiltration \
            -R {input.ref} \
            -V {params.tmpdir}/indels.vcf.gz \
            --filter-expression "QD < 2.0"              --filter-name "QD2"             \
            --filter-expression "FS > 200.0"            --filter-name "FS200"           \
            --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O {params.tmpdir}/indels_filt.vcf.gz 2>>{log}

        # --- Merge + depth/GQ filter ---
        gatk MergeVcfs \
            -I {params.tmpdir}/snps_filt.vcf.gz \
            -I {params.tmpdir}/indels_filt.vcf.gz \
            -O {params.tmpdir}/merged.vcf.gz 2>>{log}

        bcftools view \
            -f PASS \
            -e 'FORMAT/DP < {params.min_dp} || FORMAT/GQ < {params.min_gq}' \
            {params.tmpdir}/merged.vcf.gz \
        | bcftools sort \
        | bgzip -c > {output.vcf} 2>>{log}

        tabix -p vcf {output.vcf}
        rm -rf {params.tmpdir}
        """
