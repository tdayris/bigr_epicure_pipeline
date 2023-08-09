rule deeptools_alignment_sieve:
    input:
        aln="sambamba/markdup/{sample}.bam",
        aln_idx="sambamba/markdup/{sample}.bam.bai",
        blacklist=blacklist_path,
    output:
        temp("deeptools/alignment_sieve/{sample}.bam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 3,
        tmpdir=tmp,
    log:
        "logs/deeptools/alignmentsieve/{sample}.log",
    params:
        extra=lambda wildcards: get_deeptools_alignment_sieve_params(wildcards),
    wrapper:
        "master/bio/deeptools/alignmentsieve"


rule sort_deeptools_alignment_sieve:
    input:
        "deeptools/alignment_sieve/{sample}.bam",
    output:
        temp("deeptools/sorted_sieve/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{sample}.deeptools_alignment_sieve.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/sort"


rule sambamba_index_deeptools_alignment_sieve:
    input:
        "deeptools/sorted_sieve/{sample}.bam",
    output:
        temp("deeptools/sorted_sieve/{sample}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{sample}.deeptools_alignment_sieve.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/index"


rule deeptools_alignment_sieve_per_factor:
    input:
        aln="sambamba/markdup/{factor_level}.bam",
        aln_idx="sambamba/markdup/{factor_level}.bam.bai",
        blacklist=blacklist_path,
    output:
        temp("deeptools/alignment_sieve/{factor_level}.bam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 3,
        tmpdir=tmp,
    log:
        "logs/deeptools/alignmentsieve/{factor_level}.log",
    params:
        extra=lambda wildcards: get_deeptools_alignment_sieve_params(wildcards),
    wrapper:
        "master/bio/deeptools/alignmentsieve"


rule sort_deeptools_alignment_sieve_per_factor:
    input:
        "deeptools/alignment_sieve/{factor_level}.bam",
    output:
        temp("deeptools/sorted_sieve/{factor_level}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{factor_level}.deeptools_alignment_sieve.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/sort"


rule sambamba_index_deeptools_alignment_sieve_per_factor:
    input:
        "deeptools/sorted_sieve/{factor_level}.bam",
    output:
        temp("deeptools/sorted_sieve/{factor_level}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{factor_level}.deeptools_alignment_sieve.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/index"
