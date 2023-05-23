rule deeptools_compute_matrix:
    input:
        bed=expand(
            "macs2/callpeak/{sample}_peaks.{peaktype}Peak.bed",
            sample=sample_list,
            allow_missing=True
        ),
        bigwig=expand(
            "data_output/Coverage/{sample}.bw",
            sample=sample_list
        ),
        blacklist=blacklist_path,
    output:
        matrix_gz=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.gz"),
        matrix_tab=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.tab"),
        matrix_bed=temp("deeptools/matrix_files/{command}/{peaktype}.matrix.bed"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wilcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        command="{command}",
        extra="--skipZeros ",
    wrapper:
        "v1.31.1/bio/deeptools/computematrix"
