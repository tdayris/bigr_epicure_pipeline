rule link_blacklist_gz:
    input:
        config["reference"]["blacklist"],
    output:
        "reference/blacklist/{species}.{build}.{release}.bed",
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
