rule csaw_count_normalize:
    input:
        counts="csaw/filtered/{comparison_name}.RDS",
        bins="csaw/count/{comparison_name}.binned.RDS"
    output:
        rds=temp("csaw/normalized/{comparison_name}.RDS")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp
    log:
        "logs/csaw/normalize/{comparison_name}.log"
    params:
        norm_method="composition",
        extra=""
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_normalize.R"


rule csaw_count_edger:
    input:
        counts="csaw/normalized/{comparison_name}.RDS",
    output:
        csaw=temp("csaw/results/{comparison_name}.RDS"),
        edger=temp("csaw/edger/{comparison_name}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp
    log:
        "logs/csaw/edger/{comparison_name}.log"
    params:
        extra="",
        design=design.copy(),
        formula=lambda wildcards: get_edger_formula(wildcards),
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_edger.R"