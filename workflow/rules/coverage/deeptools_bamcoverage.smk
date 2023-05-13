# TODO: Make sure non-chip sonication (e.g. mnase) do not uses these parameters
# ask E. later
rule test_deeptools_bamcoverage:
    input:
        unpack(get_deeptools_bamcoverage_input),
    output:
        temp("data_output/Coverage/{sample}.bw"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        effective_genome_size=effective_genome_size,
        extra=" --normalizeUsing RPKM ",
    log:
        "logs/deeptools/coverage/{sample}.log",
    wrapper:
        "v1.30.0/bio/deeptools/bamcoverage"
