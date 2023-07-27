rule sambamba_merge_per_factors_level:
    input:
        unpack(get_sambamba_merge_per_factors_level_input_input),
    output:
        temp("sambamba/markdup/{factor_level}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/sambamba/merge/{factor_level}.log",
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/merge"


rule index_sambamba_merge_per_factors_level:
    input:
        "sambamba/markdup/{factor_level}.bam",
    output:
        "sambamba/markdup/{factor_level}.bam.bai",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/sambamba/index/{factor_level}.log",
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/index"
