rule xenome_index:
    input:
        unpack(get_xenome_index_input),
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
