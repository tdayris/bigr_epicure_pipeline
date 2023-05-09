rule macs2_callpeak_broad:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp(
            multiext(
                "macs2/callpeak/{sample}",
                "_peaks.xls",
                "_peaks.broadPeak",
                "_peaks.gappedPeak",
            )
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/macs2/callpeak/broad/{sample}.log",
    params:
        lambda wildcards: get_macs2_params(wildcards),
    wrapper:
        "v1.29.0/bio/macs2/callpeak"


rule macs2_callpeak_narrow:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp(
            multiext(
                "macs2/callpeak/{sample}",
                "_peaks.xls",
                "_peaks.narrowPeak",
                "_summits.bed",
            )
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 7,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/macs2/callpeak/narrow/{sample}.log",
    params:
        lambda wildcards: get_macs2_params(wildcards),
    wrapper:
        "v1.29.0/bio/macs2/callpeak"
