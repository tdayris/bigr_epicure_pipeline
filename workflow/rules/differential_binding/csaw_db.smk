rule csaw_count_normalize:
    input:
        counts="csaw/filtered/{model_name}.RDS",
        bins="csaw/count/{model_name}.binned.RDS",
    output:
        rds=temp("csaw/normalized/{model_name}.RDS"),
        adj_counts=temp("csaw/normalized/{model_name}.adj_counts.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/csaw/normalize/{model_name}.log",
    params:
        norm_method="composition",
        extra="",
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_normalize.R"


rule csaw_edger:
    input:
        counts="csaw/normalized/{model_name}.RDS",
        design=design_path,
    output:
        csaw=temp("csaw/results/{model_name}.RDS"),
        edger=temp("csaw/edger/{model_name}.RDS"),
        disp_png="data_output/Differential_Binding/{model_name}/Dispersion_Estimates.png",
        ql_png="data_output/Differential_Binding/{model_name}/Quasi_likelihood_dispersion.png",
        qc=temp("csaw/edger/{model_name}.sig.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/csaw/edger/{model_name}.log",
    params:
        extra="",
        formula=lambda wildcards: get_edger_formula(wildcards),
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/csaw/csaw_edger.R"
