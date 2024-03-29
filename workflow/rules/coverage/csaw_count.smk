rule csaw_fragment_length:
    input:
        design="rsamtools/design.RDS",
        read_params="csaw/readparam.{library}.RDS",
    output:
        fragment_length="csaw/fragment_length/{library}.RDS",
        png="data_output/Differential_Binding/CCF_Delay.{library}.png",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/csaw/fragment_length/{library}.log",
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_fregment_length.R"


rule csaw_window_count:
    input:
        unpack(get_csaw_count_input),
    output:
        rds=temp("csaw/count/{model_name}.{signal}.RDS"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/csaw/windowcount/{model_name}.{signal}.log",
    params:
        extra=lambda wildcards: get_csaw_count_params(wildcards),
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_count.R"


rule csaw_count_filter:
    input:
        unpack(get_csaw_filter_input),
    output:
        rds=temp("csaw/filtered/{model_name}.RDS"),
        png="data_output/Differential_Binding/{model_name}/Background_filter_abundance.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/csaw/filter/{model_name}.log",
    params:
        filter_method=lambda wildcards, input: "input_controls"
        if "input_counts" in input.keys()
        else "proportion",
        proportion=0.999,
        min_cpm=5,
        extra="",
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_filter.R"
