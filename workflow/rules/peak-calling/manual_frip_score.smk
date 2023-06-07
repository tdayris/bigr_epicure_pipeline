rule sambamba_view_compute_called_read_nb:
    input:
        "bedtools/intersect/{sample}.{peaktype}.bam",
    output:
        temp("frip/called/{sample}.{peaktype}.txt"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/frip.{sample}.{peaktype}.log",
    params:
        "--count",
    conda:
        "../../envs/sambamba.yaml"
    shell:
        "NB=$(sambamba view {params} {input}) && "
        'echo -e "{wildcards.sample}\t${{NB}}" > {output} 2> {log}'


rule sambamba_view_compute_read_nb:
    input:
        "sambamba/markdup/{sample}.bam",
    output:
        temp("frip/total/{sample}.txt"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/frip.{sample}.log",
    params:
        "--count",
    conda:
        "../../envs/sambamba.yaml"
    shell:
        "NB=$(sambamba view {params} {input}) && "
        'echo -e "{wildcards.sample}\t${{NB}}" > {output} 2> {log}'


rule manual_frip_score_concat:
    input:
        called=expand(
            "frip/called/{sample}.{peaktype}.txt",
            sample=get_tested_sample_list(),
            allow_missing=True,
        ),
        total=expand("frip/total/{sample}.txt", sample=get_tested_sample_list()),
    output:
        temp("frip/complete/{peaktype}.txt"),
    threads: 3
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/frip/complete/{peaktype}.log",
    params:
        "",
    conda:
        "../../envs/bash.yaml"
    shell:
        "paste - - <(cat {input.called}) <(cat {input.total}) > {output} 2> {log}"


rule manual_frip_score_compute:
    input:
        "frip/complete/{peaktype}.txt",
    output:
        "data_output/Peak_Calling/{peaktype}.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/compute/{peaktype}.log",
    params:
        script=workflow.source_path("../../scripts/frip/compute_frip.awk"),
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.script} {input} > {output} 2> {log}"
