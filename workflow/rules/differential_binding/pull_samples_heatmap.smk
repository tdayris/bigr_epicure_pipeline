rule deeptools_compute_matrix_per_factor:
    input:
        unpack(get_deeptools_compute_matrix_per_factor_input),
    output:
        matrix_gz=temp(
            "deeptools/matrix_files/{command}/{peaktype}.{model_name}.matrix.gz"
        ),
        matrix_tab=temp(
            "deeptools/matrix_files/{command}/{peaktype}.{model_name}.matrix.tab"
        ),
        matrix_bed=temp(
            "deeptools/matrix_files/{command}/{peaktype}.{model_name}.matrix.bed"
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
        runtime=lambda wildcards, attempt: attempt * 6 * 60,
        tmpdir=tmp,
    log:
        "logs/deeptools/computematrix/{command}.{peaktype}.{model_name}.log",
    params:
        command="{command}",
        extra=lambda wildcards: get_deeptools_compute_matrix_params(wildcards),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/computematrix"


rule deeptools_plot_heatmap_per_factor:
    input:
        "deeptools/matrix_files/{command}/{peaktype}.{model_name}.matrix.gz",
    output:
        heatmap_img=report(
            "data_output/Heatmaps/{peaktype}/{command}.{model_name}.png",
            caption="../../report/coverage/plot_heatmap.rst",
            category="Results",
            labels={"type": "figure", "category": "Results"},
        ),
        heatmap_matrix=temp(
            "deeptools/heatmap_matrix/{command}.{peaktype}.{model_name}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
        runtime=lambda wildcards, attempt: attempt * 6,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_heatmap/{peaktype}.{command}.{model_name}.log",
    params:
        extra=lambda wildcards: get_deeptools_plotheatmap_params(wildcards),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotheatmap"
