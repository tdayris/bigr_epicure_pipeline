rule deeptools_plot_fingerprint:
    input:
        unpack(get_deeptools_plotfingerprint_input),
    output:
        fingerprint=report(
            "data_output/QC/fingerprint.png",
            caption="../../report/coverage/plot_fingerprint.rst",
            category="Quality Control",
            labels={
                "type": "figure",
                "category": "QC",
            },
        ),
        counts=temp("deeptools/plot_fingerprint/raw_counts.tab"),
        qc_metrics=temp("deeptools/plot_fingerprint/qc_metrics.txt"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_fingerprint.log",
    params:
        lambda wildcards: get_deeptools_fingerprint_params(wildcards),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotfingerprint"
