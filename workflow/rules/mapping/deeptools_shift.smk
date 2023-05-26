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
        "logs/deeptools/alignmentsieve/{sample}.log",
    params:
        extra=" --ATACshift ",
    wrapper:
        "master/bio/deeptools/alignmentsieve"
