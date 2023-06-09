rule sambamba_sort_filtered:
    input:
        unpack(get_sambamba_sort_filtered),
    output:
        temp("sambamba/filtered/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{sample}.bowtie2cr.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/sort"


rule sambamba_markdup:
    input:
        "sambamba/filtered/{sample}.bam",
    output:
        temp("sambamba/markdup/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra=lambda wildcards: get_sambamba_markdup_params(wildcards),
    log:
        "logs/sambamba/markdup/{sample}.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/markdup"
