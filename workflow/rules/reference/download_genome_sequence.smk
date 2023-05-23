rule get_genome:
    output:
        "reference/{species}.{build}.{release}.fasta",
    params:
        species="{species}",
        datatype="dna",
        build="{release}",
        release="{build}",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/reference/sequence/{species}.{build}.{release}.log",
    cache: "omit-software"
    wrapper:
        "v1.31.1/bio/reference/ensembl-sequence"
