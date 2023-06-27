rule sambamba_merge_per_factor:
    input:
        unpack(get_sambamba_merge_input),
    output:
        temp("sambamba/merged/{factor_level}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 6,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/sambamba/merge/{factor_level}.log",
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/merge"


rule sambamba_index_per_factor:
    input:
        "sambamba/merged/{factor_level}.bam",
    output:
        temp("sambamba/merged/{factor_level}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{factor_level}.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/index"


rule deeptools_bamcoverage_per_factor:
    input:
        bam="sambamba/merged/{factor_level}.bam",
        bai="sambamba/merged/{factor_level}.bam.bai",
        blacklist=blacklist_path,
    output:
        "data_output/Coverage/{factor_level}.bw",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        effective_genome_size=effective_genome_size,
        extra=lambda wildcards: get_deeptools_bamcoverage_params(wildcards),
    log:
        "logs/deeptools/coverage/{factor_level}.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/bamcoverage"
