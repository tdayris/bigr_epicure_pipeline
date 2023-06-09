rule medips_edger_diff_peak:
    input:
        rds="medips/meth/{model_name}.RDS",
    output:
        rds="medips/edger/{model_name}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/edger/{model_name}.log",
    params:
        extra="",
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medip-seq/medips_edger.R"
