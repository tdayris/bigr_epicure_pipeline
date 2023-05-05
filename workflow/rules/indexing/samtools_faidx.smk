rule samtools_index:
    input:
        "reference/{species}.{build}.{release}.fasta",
    output:
        "reference/{species}.{build}.{release}.fasta.fai",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1 * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/samtools/faidx/{species}.{build}.{release}.log",
    params:
        extra="",
    wrapper:
        "v1.29.0/bio/samtools/faidx"
