rule medips_import_sample_bam:
    input:
        unpack(get_medips_import_sample_bam_input),
    output:
        rds=temp("medips/import/{comparison_name}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/import/{comparison_name}.log",
    params:
        extra=lambda wildcards: get_medips_params_extra(wildcards),
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medips_load.R"


rule medips_make_coupling_vector:
    input:
        rds="medips/import/{comparison_name}.RDS",
    output:
        rds="medips/coupling/{comparison_name}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/coupling/{comparison_name}.log",
    params:
        pattern="CG",
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medips_coupling.R"
