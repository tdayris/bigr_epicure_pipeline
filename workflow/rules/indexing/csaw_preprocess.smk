rule csaw_readparam:
    input:
        blacklist=blacklist_path,
    output:
        rds=temp("csaw/readparam.{library}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/csaw/readparam/{library}.log",
    params:
        extra=lambda wildcards: get_csaw_read_param(wildcards),
        organism=build,
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_readparam.R"
