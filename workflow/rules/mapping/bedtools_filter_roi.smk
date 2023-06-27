rule bedtools_filter_roi:
    input:
        unpack(get_bedtools_filter_roi_input),
    output:
        temp("bedtools/filter_roi/{sample}.bam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/bedtools/filter_roi/{sample}.log",
    params:
        extra=" -wa -sorted ",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "bedtools intersect "
        "{params.extra} "
        "-abam {input.left} "
        "-b {input.right} "
        "> {output} 2> {log}"


rule sambamba_index_roi:
    input:
        "bedtools/filter_roi/{sample}.bam",
    output:
        temp("bedtools/filter_roi/{sample}.bam.bai"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="",
    log:
        "logs/sambamba/index/{sample}.raw.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/index"
