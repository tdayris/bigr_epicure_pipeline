rule medips_meth_coverage:
    input:
        unpack(get_medips_meth_coverage_input),
    output:
        rds=temp("medips/meth/{model_name}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/meth/{model_name}.log",
    params:
        extra="",
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medip-seq/medips_meth.R"
