rule link_bowtie2_index:
    input:
        config["reference"]["bowtie2_index"],
    output:
        bowtie2_index_path,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/ln/blacklist/{species}.{build}.{release}.log",
    params:
        extra_ln="--relative --verbose --symbolic --force",
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/misc/rename.py"
