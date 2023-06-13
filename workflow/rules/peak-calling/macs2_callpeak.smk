rule macs2_callpeak_broad:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp("macs2/callpeak_broad/{sample}_peaks.xls"),
        temp(ensure("macs2/callpeak_broad/{sample}_peaks.broadPeak", non_empty=True)),
        temp("macs2/callpeak_broad/{sample}_peaks.gappedPeak"),
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
        "v1.32.1/bio/macs2/callpeak"


rule macs2_save_broad:
    input:
        "macs2/callpeak_broad/{sample}_peaks.xls",
    output:
        "data_output/Peak_Calling/broad/Macs2/{sample}_peaks.xls",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/rsync/macs2/{sample}.broad.log",
    conda:
        "../../envs/bash.yaml"
    params:
        "-cvhP",
    shell:
        "rsync {params} {input} {output} > {log} 2>&1"


rule macs2_callpeak_narrow:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp("macs2/callpeak_narrow/{sample}_peaks.xls"),
        temp(ensure("macs2/callpeak_narrow/{sample}_peaks.narrowPeak", non_empty=True)),
        temp("macs2/callpeak_narrow/{sample}_summits.bed"),
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
        "v1.32.1/bio/macs2/callpeak"


rule macs2_save_narrow:
    input:
        "macs2/callpeak_narrow/{sample}_peaks.xls",
    output:
        "data_output/Peak_Calling/narrow/Macs2/{sample}_peaks.xls",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/rsync/macs2/{sample}.narrow.log",
    conda:
        "../../envs/bash.yaml"
    params:
        "-cvhP",
    shell:
        "rsync {params} {input} {output} > {log} 2>&1"


# macs2_mode = {narrow, broad, gapped}
rule macs2_peaks_to_bed:
    input:
        peak="macs2/callpeak_{macs2_mode}/{sample}_peaks.{macs2_mode}Peak",
    output:
        bed=temp("macs2/callpeak_{macs2_mode}/{sample}_peaks.{macs2_mode}Peak.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/macs2/peak2bed/{sample}.{macs2_mode}.log",
    params:
        script=workflow.source_path("../../scripts/formatter/macs2_to_bed6.awk"),
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.script} {input} > {output} 2> {log}"
