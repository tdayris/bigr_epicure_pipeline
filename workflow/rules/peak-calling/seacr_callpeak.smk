rule bedtools_genomecoveragebed:
    input:
        "sambamba/markdup/{sample}.bam",
    output:
        temp("bedtools/genomecov/{sample}.bg"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        "logs/bedtools/genomecov/{sample}.log",
    params:
        extra=" -bg ",
    wrapper:
        "v1.32.1/bio/bedtools/genomecov"


rule seacr_callpeak:
    input:
        unpack(get_seacr_callpeak_input),
    output:
        temp("seacr/raw/{sample}.{seacr_mode}.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 60 * 1,
        tmpdir=tmp,
    log:
        "logs/seacr/{sample}.{seacr_mode}.log",
    params:
        extra=lambda wildcards: get_seacr_callpeak_params(wildcards),
    conda:
        "../../envs/seacr.yaml"
    shell:
        "SEACR_1.3.sh "
        "{input.exp_bg} "
        "{params.extra} "
        "> {log} 2>&1"


rule format_seacr_output:
    input:
        "seacr/raw/{sample}.{seacr_mode}.bed",
    output:
        temp("seacr/{seacr_mode}/{sample}.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/seacr/format_awk/{sample}.{seacr_mode}.log",
    params:
        script=workflow.source_path("../../scripts/formatter/seacr_to_bed6.awk"),
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.script} {input} > {output} 2> {log}"
