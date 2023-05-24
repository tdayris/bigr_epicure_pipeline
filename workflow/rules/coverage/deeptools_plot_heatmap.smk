rule deeptools_plot_heatmap:
    input:
        "deeptools/matrix_files/{command}/{peaktype}.matrix.gz",
    output:
        report(
            "data_output/Heatmaps/{peaktype}/{command}.png",
            caption="../../report/coverage/plot_heatmap.rst",
            category="Results",
            labels={
                "type": "figure",
                "category": "Results"
            }
        )
    threads: 1
    resources:
        mem_mb=lambda wilcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_heatmap/{peaktype}.{command}.log",
    params:
        extra="--plotFileFormat 'png'",
    wrapper:
        "v1.31.1/bio/deeptools/plotheatmap"
