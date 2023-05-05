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
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools.log",
    params:
        extra="--ATACshift",
    wrapper:
        "master/bio/deeptools/alignmentsieve"


rule sambamba_index_deeptools_alignment_sieve:
    input:
        "deeptools/alignment_sieve/{sample}.bam",
    output:
        temp("deeptools/alignment_sieve/{sample}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{sample}.shifted.log",
    wrapper:
        "v1.29.0/bio/sambamba/index"
