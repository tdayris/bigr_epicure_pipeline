rule deeptools_bampe_fragment_size:
    input:
        bam=expand(
            "sambamba/filtered/{sample}.bam",
            sample=get_tested_sample_list(),
        ),
        blacklist=blacklist_path,
    output:
        hist="data_output/QC/Fragment_Size/{sample}.png",
        metrics=temp("deeptools/bampe_fs/{sample}_metrics.txt"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        "logs/deeptools/bampe_fs/{sample}.log",
    params:
        extra=lambda wildcards: get_deeptools_bampefragmentsize_params(wildcards),
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "bamPEFragmentSize "
        "{params.extra} "
        "--bamfiles {input.bam} "
        "--histogram {output.hist} "
        "--numberOfProcessors {threads} "
        "--blackListFileName {input.blacklist} "
        "> {output.metrics} 2> {log} "
