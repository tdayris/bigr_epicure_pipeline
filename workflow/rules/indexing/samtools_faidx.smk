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
        f"{snakemake_wrappers_version}/bio/samtools/faidx"


rule picard_create_sequence_dict:
    input:
        "reference/{species}.{build}.{release}.fasta",
    output:
        "reference/{species}.{build}.{release}.dict",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/picard/create_dict/{species}.{build}.{release}.log",
    params:
        extra="",
    wrapper:
        f"{snakemake_wrappers_version}/bio/picard/createsequencedictionary"


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
        extra="",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/index"
