rule rename_input_files_single:
    output:
        temp("data_input/{sample}.fq.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/rename_concat/{sample}.log",
    params:
        input=lambda wildcards: get_rename_input(wildcards, design),
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/misc/rename.py"


rule rename_input_files_paired:
    output:
        temp("data_input/{sample}.{stream}.fq.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    message:
        "output: {output}, params: {params.input}, wildcards: {wildcards.sample}, {wildcards.stream}"
    log:
        "logs/rename_concat/{sample}.{stream}.log",
    params:
        input=lambda wildcards: get_rename_input(wildcards, design),
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/misc/rename.py"
