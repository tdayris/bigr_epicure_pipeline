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


rule macs2_save_broad:
    input:
        "macs2/callpeak/{sample}_peaks.xls",
    output:
        "data_output/Peak_Calling/macs2/{sample}_broad_peaks.xls",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 128,
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


rule macs2_save_narrow:
    input:
        "macs2/callpeak/{sample}_peaks.xls",
    output:
        "data_output/Peak_Calling/macs2/{sample}_narrow_peaks.xls",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 128,
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


# peaktype = {narrow, broad, gapped}
rule macs2_peaks_to_bed:
    input:
        peak="macs2/callpeak/{sample}_peaks.{peaktype}Peak",
    output:
        bed=temp("macs2/callpeak/{sample}_peaks.{peaktype}Peak.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/macs2/peak2bed/{sample}.{peaktype}.log",
    params:
        script=workflow.source_path("../../scripts/formatter/macs2_to_bed6.awk"),
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.script} {input} > {output} 2> {log}"
