rule picard_collect_multiple_metrics:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bai="sambamba/markdup/{sample}.bam.bai",
        ref=genome_fasta_path,
        ref_fai=genome_fai_path,
        ref_dict=genome_dict_path,
    output:
        temp(
            multiext(
                "picard/collectmultiplemetrics/stats/{sample}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            )
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/picard/multiple_metrics/{sample}.log",
    params:
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    wrapper:
        "v1.31.1/bio/picard/collectmultiplemetrics"
