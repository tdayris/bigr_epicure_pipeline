rule fasta_to_two_bit:
    input:
        genome_fasta_path,
    output:
        genome_twobit_path,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/ucsc/fatotwobit.log",
    params:
        extra="",
    wrapper:
        "v1.32.1/bio/ucsc/faToTwoBit"
