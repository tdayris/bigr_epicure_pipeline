rule deeptools_plot_correlation:
    input:
        bw="deeptools/multibigwig/summary.bw",
    output:
        png="data_output/QC/Correlation.{plot_type}.png",
        stats=temp("deeptools/plot_corr/{plot_type}.stats"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_corr/summary/{plot_type}.log",
    params:
        extra=lambda wildcards: get_deeptools_plot_correlation_params(wildcards),
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "plotCorrelation --corData {input.bw} "
        "--plotFile {output.png} "
        "--outFileCorMatrix {output.stats} "
        "{params.extra} > {log} 2>&1 "
