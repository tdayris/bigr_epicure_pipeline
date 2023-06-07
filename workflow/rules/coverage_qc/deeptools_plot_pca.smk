rule deeptools_plot_pca:
    input:
        bw="deeptools/multibigwig/summary.bw",
    output:
        png="data_output/QC/PCA.png",
        stats=temp("deeptools/plot_pca.stats"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/deeptools/plot_pca/summary.log",
    params:
        "--plotFileFormat png ",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "plotPCA {params} --corData {input.bw} "
        "--plotFile {output.png} "
        "--outFileNameData {output.stats} "
        "> {log} 2>&1 "
