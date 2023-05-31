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
        extra=lambda wilcards: get_csaw_read_param(wilcards)
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_readparam.R"


rule csaw_window_count:
    input:
        unpack(get_csaw_count_input)
    output:
        rds=temp("csaw/count/{comparison_name}.{signal}.RDS")
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp
    log:
        "logs/csaw/windowcount/{comparison_name}.{signal}.log"
    params:
        extra=lambda wilcards: get_csaw_count_params(wilcards)
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_count.R"


rule csaw_count_filter:
    input:
        unpack(get_csaw_filter_input)
    output:
        rds=temp("csaw/filtered/{comparison_name}.RDS")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp
    log:
        "logs/csaw/filter/{comparison_name}.log"
    params:
        filter_method=lambda wildcards, input: "input_controls" if "input_counts" in input.keys() else "proportion",
        proportion=0.999,
        min_cpm=5,
        extra=""
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_filter.R"