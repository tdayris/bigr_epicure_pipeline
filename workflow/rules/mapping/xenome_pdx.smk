rule xenome_classify_pairs:
    input:
        unpack(get_xenome_classify_input),
    output:
        temp(
            expand(
                "xenome/classify/{sample}_{origin}_{stream}.fastq",
                origin=["ambiguous", "both", "graft", "host", "neither"],
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 15,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        general="logs/xenome/classify/{sample}.pairs.log",
        tool="logs/xenome/classify/{sample}.pairs.internal.log",
    params:
        extra=" --verbose --pairs ",
        out_prefix=lambda wildcards, output: output[0][: -len("_ambiguous_1.fastq")],
        idx_prefix=lambda wildcards, input: sorted(input.index)[0][
            : -len("-both.header")
        ],
    envmodules:
        "xenome/1.0.0_patched",
    shell:
        "xenome classify {params.extra} "
        "--log-file {log.tool} "
        "--num-threads {threads} "
        "--tmp-dir {resources.tmpdir} "
        "--max-memory {resources.mem_mb} "
        "--output-filename-prefix {params.out_prefix} "
        "--prefix {params.idx_prefix} "
        "--fastq-in {input.r1} --fastq-in {input.r2} "
        "> {log.general} 2>&1 "


rule xenome_classify_single:
    input:
        unpack(get_xenome_classify_input),
    output:
        temp(
            expand(
                "xenome/classify/{sample}_{origin}.fastq",
                origin=["ambiguous", "both", "graft", "host", "neither"],
                allow_missing=True,
            )
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 15,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        general="logs/xenome/classify/{sample}.pairs.log",
        tool="logs/xenome/classify/{sample}.pairs.internal.log",
    params:
        extra=" --verbose --pairs ",
        out_prefix=lambda wildcards, output: output[0][: -len("_ambiguous.fastq")],
        idx_prefix=lambda wildcards, input: sorted(input.index)[0][
            : -len("-both.header")
        ],
    envmodules:
        "xenome/1.0.0_patched",
    shell:
        "xenome classify {params.extra} "
        "--log-file {log.tool} "
        "--num-threads {threads} "
        "--tmp-dir {resources.tmpdir} "
        "--max-memory {resources.mem_mb} "
        "--output-filename-prefix {params.out_prefix} "
        "--prefix {params.idx_prefix} "
        "--fastq-in {input.reads} "
        "> {log.general} 2>&1 "
