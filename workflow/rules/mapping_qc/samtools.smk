rule samtools_stats:
    input:
        unpack(get_samtools_stats_input)
    output:
        "samtools/stats/{sample}.{step}.txt",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
        region="",
    log:
        "logs/samtools/stats/{sample}.{step}.log",
    wrapper:
        "v1.29.0/bio/samtools/stats"