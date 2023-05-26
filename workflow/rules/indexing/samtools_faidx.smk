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
        "v1.31.1/bio/samtools/faidx"


rule sambamba_bam_index:
    input:
        "{tool}/{subcommand}/{sample}.bam",
    output:
        "{tool}/{subcommand}/{sample}.bam.bai",
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/samtools/index/{tool}.{subcommand}/{sample}.log",
    params:
        "",
    wrapper:
        "v1.31.1/bio/sambamba/index"
