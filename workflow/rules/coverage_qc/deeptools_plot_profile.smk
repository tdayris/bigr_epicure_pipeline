rule deeptools_plot_profile:
    input:
        "deeptools/matrix_files/{command}/{peaktype}.matrix.gz",
    output:
        plot_img="data_output/Peak_Calling/{peaktype}/{command}.png",
        data=temp("deeptools/plot_profile/{peaktype}/{command}.tab"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_pca/summary/{peaktype}.{command}.log",
    params:
        lambda wildcards: get_deeptools_bampefragmentsize_params(wildcards),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotprofile"
