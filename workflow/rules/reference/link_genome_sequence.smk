rule link_genome_sequence:
    input:
        config["reference"]["genome_fasta"],
    output:
        "reference/{species}.{build}.{release}.fasta",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/bigr_copy/blacklist/{species}.{build}.{release}.log",
    params:
        extra_ln="--relative --verbose --symbolic --force",
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/misc/rename.py"
