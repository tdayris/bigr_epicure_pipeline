rule sambamba_sort_filtered:
    input:
        "samtools/view/{sample}.bam",
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
        "v1.31.1/bio/sambamba/sort"


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
        extra="--remove-duplicates --overflow-list-size 600000",
    log:
        "logs/sambamba/markdup/{sample}.log",
    wrapper:
        "v1.31.1/bio/sambamba/markdup"
