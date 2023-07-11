rule rename_input_files_single:
    input:
        unpack(get_rename_input),
    output:
        temp("data_input/reads/{sample}.fastq.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/rename_concat/{sample}.log",
    params:
        "",
    conda:
        "../../envs/python.yaml"
    script:
        "../../script/misc/rename.py"


rule rename_input_files_paired:
    input:
        unpack(get_rename_input),
    output:
        temp(
            expand(
                "data_input/reads/{sample}.{stream}.fastq.gz",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/rename_concat/{sample}.log",
    params:
        "",
    conda:
        "../../envs/python.yaml"
    script:
        "../../script/misc/rename.py"
