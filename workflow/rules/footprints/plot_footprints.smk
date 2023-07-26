rule plot_footprints:
    input:
        unpack(get_plot_footprints_input),
    output:
        png="data_output/Motifs/Fingerprints/{sample}/{sample}.{motif}.png",
        bam=temp("fingerprints/{sample}.{motif}.bam"),
        bai=temp("fingerprints/{sample}.{motif}.bam.bai"),
        rda=temp("fingerprints/{sample}.{motif}.RData"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        "logs/plot_footprints/{sample}.log"
    params:
        name="{sample}",
        motif="{motif}",
    conda:
        "../../envs/factorfootprints.yaml"
    script:
        "../../scripts/factorfootprints/factor_footprints.R"
