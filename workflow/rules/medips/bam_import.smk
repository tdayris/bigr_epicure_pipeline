rule medips_import_sample_bam:
    input:
        bam=expand(
            "sambamba/markdup/{sample}.bam",
            sample=get_samples_per_condition(wildcards),
        ),
        bai=expand(
            "sambamba/markdup/{sample}.bam.bai",
            sample=get_samples_per_condition(wildcards),
        ),
    output:
        rds=temp("medips/import/{condition}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/import/{condition}.log",
    params:
        extra=lambda wildcards: get_medips_params_extra(wildcards),
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medips_load.R"


rule medips_make_coupling_vector:
    input:
        rds=temp("medips/import/{condition}.RDS"),
    output:
        rds=temp("medips/coupling/{condition}.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/medips/coupling/{condition}.log",
    params:
        pattern="CG",
    conda:
        "../../envs/medips.yaml"
    script:
        "../../scripts/medips_coupling.R"
