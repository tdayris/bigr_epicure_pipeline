rule bowtie2_align:
    input:
        unpack(get_bowtie2_align_input),
    output:
        temp("bowtie2/align/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/bowtie2/align/{sample}.log",
    params:
        extra="--very-sensitive",
    wrapper:
        "v1.32.1/bio/bowtie2/align"


rule sambamba_sort_bowtie2_aligned:
    input:
        "bowtie2/align/{sample}.bam",
    output:
        temp("sambamba/sort/{sample}.raw.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{sample}.bowtie2cr.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        "v1.31.1/bio/sambamba/sort"


rule sambamba_index_raw_bowtie2:
    input:
        "sambamba/sort/{sample}.raw.bam",
    output:
        temp("sambamba/sort/{sample}.raw.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{sample}.raw.log",
    wrapper:
        "v1.31.1/bio/sambamba/index"
