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
        "v1.29.0/bio/bowtie2/align"


rule sambamba_sort_bowtie2_aligned:
    input:
        "bowtie2/align/{sample}.bam",
    output:
        temp("sambamba/sort/{sample}.raw.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{sample}.bowtie2.log",
    params:
        extra=lambda wildcards, attempt: f"--memory-limit {attempt * 19 * 1024}MiB",
    wrapper:
        "1.29.0/bio/sambamba/sort"


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
        "v1.29.0/bio/sambamba/index"
