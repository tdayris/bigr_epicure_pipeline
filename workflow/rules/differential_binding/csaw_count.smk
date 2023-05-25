rule csaw_window_count:
    input:
        unpack(get_csaw_count_input)
    output:
        rds=temp("csaw/count/{sample}.RDS")
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp
    log:
        "logs/csaw/windowcount/{sample}.log"
    params:
        read_param=lambda wildcards: get_read_params()