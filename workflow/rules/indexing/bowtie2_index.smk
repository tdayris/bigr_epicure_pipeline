rule bowtie2_build:
    input:
        ref="reference/{species}.{build}.{release}.fasta",
    output:
        multiext(
            "reference/bowtie2_index/{species}.{build}.{release}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    cache: True
    log:
        "logs/bowtie2/build/{species}.{build}.{release}.log",
    params:
        extra="-f",
    wrapper:
        "v1.32.1/bio/bowtie2/build"
