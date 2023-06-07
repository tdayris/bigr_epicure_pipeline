rule fastp_single_ended:
    input:
        unpack(get_fastp_input),
    output:
        trimmed=temp("fastp/trimmed/se/{sample}.fastq"),
        html="data_output/QC/fastp/{sample}.se.html",
        json=temp("fastp/report/se/{sample}.fastp.json"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/fastp/{sample}.se.log",
    params:
        adapters="",
        extra=lambda wildcards: get_fastp_params(wildcards),
    wrapper:
        "v1.31.1/bio/fastp"


rule fastp_pair_ended:
    input:
        unpack(get_fastp_input),
    output:
        trimmed=temp(
            expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html="data_output/QC/fastp/{sample}.pe.html",
        json=temp("fastp/report/pe/{sample}.fastp.json"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/fastp/{sample}.pe.log",
    params:
        adapters="",
        extra=lambda wildcards: get_fastp_params(wildcards),
    wrapper:
        "v1.31.1/bio/fastp"
