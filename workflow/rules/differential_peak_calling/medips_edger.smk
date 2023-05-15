rule medips_edger_diff_peak:
    input:
        rds="medips/meth/{comparison}.RDS",
    output:
        rds="medips/edger/{comparison}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/edger/{comparison}.log",
    params:
        extra="",
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medips_edger.R"
