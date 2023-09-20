rule link_xenome_index:
    input:
        config["reference"]["xenome_index"],
    output:
        multiext(
        "reference/xenome/index/pdx-both",
            ".idx-both.header",
            ".idx-both.kmers-d0",
            ".idx-both.kmers-d1",
            ".idx-both.kmers.header",
            ".idx-both.kmers.high-bits",
            ".idx-both.kmers.low-bits.lwr",
            ".idx-both.kmers.low-bits.upr",
            ".idx-both.lhs-bits",
            ".idx-both.rhs-bits",
            ".idx-graft.header",
            ".idx-graft.kmers-d0",
            ".idx-graft.kmers-d1",
            ".idx-graft.kmers.header",
            ".idx-graft.kmers.high-bits",
            ".idx-graft.kmers.low-bits.lwr",
            ".idx-graft.kmers.low-bits.upr",
            ".idx-host.header",
            ".idx-host.kmers-d0",
            ".idx-host.kmers-d1",
            ".idx-host.kmers.header",
            ".idx-host.kmers.high-bits",
            ".idx-host.kmers.low-bits.lwr",
            ".idx-host.kmers.low-bits.upr",
        )
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/bigr_copy/xenome/{species}.{build}.{release}.log",
    params:
        extra_ln="--relative --verbose --symbolic --force",
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/misc/rename.py"