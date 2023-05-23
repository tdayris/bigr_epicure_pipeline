rule download_genome_annotation:
    output:
        "reference/{species}.{build}.{release}.gtf",
    params:
        species="{species}",
        release="{release}",
        build="{build}",
        flavor="",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/reference/annotation/{species}.{build}.{release}.log",
    cache: "omit-software"
    wrapper:
        "v1.31.1/bio/reference/ensembl-annotation"
