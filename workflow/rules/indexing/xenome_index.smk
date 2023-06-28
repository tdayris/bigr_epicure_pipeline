rule xenome_index:
    input:
        unpack(get_xenome_index_input),
    output:
        multiext(
            "xenome/index/pdx-both",
            ".header",
            ".kmers-d0",
            ".kmers-d1",
            ".kmers.header",
            ".kmers.high-bits",
            ".kmers.low-bits.lwr",
            ".kmers.low-bits.upr",
            ".lhs-bits",
            ".rhs-bits",
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 15,
        runtime=lambda wildcards, attempt: attempt * 60 * 12,
        tmpdir=tmp,
    cache: True
    log:
        general="logs/xenome/index/main.log",
        tool="logs/xenome/index/command.log",
    params:
        extra=" --verbose ",
        prefix=lambda wildcards, output: output[0].split("-")[0],
    envmodules:
        "xenome/1.0.0_patched",
    shell:
        "xenome index {params.extra} "
        "--num-threads {threads} "
        "--max-memory {resources.mem_mb} "
        "--prefix {params.prefix} "
        "--graft {input.human} "
        "--host {input.mouse} "
        "--log-file {log.tool} "
        "--tmp-dir {resources.tmpdir} "
        "> {log.general} 2>&1 "
