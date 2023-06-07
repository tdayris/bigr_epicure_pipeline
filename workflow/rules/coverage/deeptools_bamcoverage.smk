# TODO: Make sure non-chip sonication (e.g. mnase) do not uses these parameters
# ask E. later
rule deeptools_bamcoverage:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bai="sambamba/markdup/{sample}.bam.bai",
        blacklist=blacklist_path,
    output:
        "data_output/Coverage/{sample}.bw",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        effective_genome_size=effective_genome_size,
        extra=" --normalizeUsing RPKM --binSize 50",
    log:
        "logs/deeptools/coverage/{sample}.log",
    wrapper:
        "v1.31.1/bio/deeptools/bamcoverage"


rule deeptools_multibigwig_summary:
    input:
        bw=expand("data_output/Coverage/{sample}.bw", sample=get_tested_sample_list()),
        blacklist=blacklist_path,
    output:
        bw=temp("deeptools/multibigwig/summary.bw"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 12,
        runtime=lambda wildcards, attempt: attempt * 60 * 3,
        tmpdir=tmp,
    log:
        "logs/deeptools/multibigwig/summary.log",
    params:
        "",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins {params} "
        "--bwfiles {input.bw} "
        "--outFileName {output.bw} "
        "--blackListFileName {input.blacklist} "
        "--numberOfProcessors {threads} "
        "> {log} 2>&1 "


rule deeptools_plotcoverage:
    input:
        bam=expand("sambamba/markdup/{sample}.bam", sample=design.index),
        bai=expand("sambamba/markdup/{sample}.bam.bai", sample=design.index),
        blacklist=blacklist_path,
    output:
        plot=report(
            "data_output/Coverage/PlotCoverage.png",
            caption="../../report/coverage/plot_coverage.rst",
            category="Quality Control",
            labels={
                "type": "figure",
                "category": "QC",
            },
        ),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/plotcoverage.log",
    params:
        extra="--skipZeros --centerReads --ignoreDuplicates --minMappingQuality 10",
    wrapper:
        "v1.31.1/bio/deeptools/plotcoverage"
