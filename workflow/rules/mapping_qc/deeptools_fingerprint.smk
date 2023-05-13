rule deeptools_plot_fingerprint:
    input:
        unpack(get_deeptools_plotfingerprint_input)
    output:
        fingerprint="data_output/QC/fingerprint.png",
        counts=temp("deeptools/plot_fingerprint/raw_counts.tab"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_fingerprint.log"
    params:
        " --skipZeros "
    wrapper:
        "v1.29.0/bio/deeptools/plotfingerprint"