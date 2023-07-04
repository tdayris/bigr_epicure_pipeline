rule rsamtools_bam_file:
    input:
        unpack(get_rsamtools_bam_file),
    output:
        rds="rsamtools/{sample}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/rsamtools/{sample}.log",
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/rsamtools/create_bam_file.R"


rule rsamtools_design:
    input:
        bam_files=expand("rsamtools/{sample}.RDS", sample=get_tested_sample_list()),
        design=design_path,
    output:
        qc=temp("rsamtools/qc.tsv"),
        design=temp("rsamtools/design.RDS"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/rsamtools/design.log",
    conda:
        "../../envs/csaw.yaml"
    script:
        "../../scripts/rsamtools/create_design.R"
