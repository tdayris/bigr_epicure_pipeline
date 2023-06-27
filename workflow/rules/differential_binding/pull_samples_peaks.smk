rule cat_bed_per_factor:
    input:
        unpack(get_bedtools_merge_factor_input),
    output:
        temp("bigr/cat_bed/{factor_level}.{peaktype}.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    group:
        "bed_per_factor"
    log:
        "logs/bigr/cat/{factor_level}.{peaktype}.log",
    params:
        "",
    conda:
        "../../envs/bash.yaml"
    shell:
        "cat {params} {input} > {output} 2> {log}"


rule sort_bed_per_factor:
    input:
        "bigr/cat_bed/{factor_level}.{peaktype}.bed",
    output:
        temp("bigr/cat_bed/{factor_level}.{peaktype}.sorted.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    group:
        "bed_per_factor"
    log:
        "logs/bigr/sort/{factor_level}.{peaktype}.log",
    params:
        "-k1,1 -k2,2n",
    conda:
        "../../envs/bash.yaml"
    shell:
        "sort {params} {input} > {output} 2> {log}"


rule bedtools_merge_per_factor:
    input:
        "bigr/cat_bed/{factor_level}.{peaktype}.sorted.bed",
    output:
        temp("bedtools/merge/{factor_level}.{peaktype}.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/bedtools/merge/{factor_level}.{peaktype}.log",
    params:
        extra=" -d 5 ",
    wrapper:
        "v2.0.0/bio/bedtools/merge"
