rule deeptools_compute_matrix:
    input:
        unpack(get_deeptools_compute_matrix_input),
    output:
        matrix_gz=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.gz"),
        matrix_tab=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.tab"),
        matrix_bed=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.bed"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
        runtime=lambda wildcards, attempt: attempt * 6 * 60,
        tmpdir=tmp,
    log:
        "logs/deeptools/computematrix/{command}.{peaktype}.log",
    params:
        command="{command}",
        extra=lambda wildcards: get_deeptools_compute_matrix_params(wildcards),
    wrapper:
        "v1.32.1/bio/deeptools/computematrix"
