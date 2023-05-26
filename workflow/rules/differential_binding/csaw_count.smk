rule csaw_readparam:
    output:
        rds=temp("csaw/readparam.{library}.RDS")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp
    log:
        "logs/csaw/readparam.log"
    params:
        extra=lambda wilcards: "minq=30, pe='both'" if str(wilcards.library) == "pe" else "minq=30"
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_readparam.R"


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
        extra=""
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_count.R"