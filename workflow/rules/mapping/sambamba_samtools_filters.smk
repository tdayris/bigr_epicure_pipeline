rule sambamba_quality_filter:
    input:
        "sambamba/sort/{sample}.raw.bam" "sambamba/sort/{sample}.raw.bam.bai",
    output:
        temp("sambamba/view/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/{sample}.filter.log",
    params:
        extra="--with-header --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)'",
    wrapper:
        "1.29.0/bio/sambamba/view"


rule samtools_filter_canonical_chromosomes:
    input:
        "sambamba/view/{sample}.bam",
    output:
        bam=temp("samtools/view/{sample}.bam"),
        idx=temp("samtools/view/{sample}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "{sample}.log",
    params:
        extra="",
        region=" ".join(canonical_chromosomes),
    wrapper:
        "v1.29.0/bio/samtools/view"
