rule deeptools_plot_heatmap:
    input:
        "deeptools/matrix_files/{command}/{peaktype}.matrix.gz",
    output:
        heatmap_img=report(
            "data_output/Heatmaps/{peaktype}/{command}.png",
            caption="../../report/coverage/plot_heatmap.rst",
            category="Results",
            labels={"type": "figure", "category": "Results"},
        ),
        heatmap_matrix=temp("deeptools/heatmap_matrix/{command}.{peaktype}.tab"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
        runtime=lambda wildcards, attempt: attempt * 6,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_heatmap/{peaktype}.{command}.log",
    params:
        extra="--plotFileFormat 'png'",
    wrapper:
        "v1.32.1/bio/deeptools/plotheatmap"
