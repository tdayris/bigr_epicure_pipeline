# TODO: Make sure non-chip sonication (e.g. mnase) do not uses these parameters
# ask E. later
rule test_deeptools_bamcoverage:
    input:
        unpack(get_deeptools_bamcoverage_input),
    output:
        temp("data_output/Coverage/{sample}.bw"),
    params:
        effective_genome_size=1000,
        extra="--normalizeUsing RPKM",
    log:
        "logs/coverage.log",
    wrapper:
        "master/bio/deeptools/bamcoverage"
