# TODO: Make sure non-chip sonication (e.g. mnase) do not uses these parameters
# ask E. later
rule deeptools_bamcoverage:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bai="sambamba/markdup/{sample}.bam.bai",
        blacklist=blacklist_path,
    output:
        temp("data_output/Coverage/{sample}.bw"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        effective_genome_size=effective_genome_size,
        extra=" --normalizeUsing RPKM ",
    log:
        "logs/deeptools/coverage/{sample}.log",
    wrapper:
        "v1.31.1/bio/deeptools/bamcoverage"


rule deeptools_plotcoverage:
    input:
        bam=expand("sambamba/markdup/{sample}.bam", sample=sample_list),
        bai=expand("sambamba/markdup/{sample}.bam.bai", sample=sample_list),
        blacklist=blacklist_path,
    output:
        plot="data_output/Coverage/PlotCoverage.png"
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/plotcoverage.log"
    params:
        extra="--skipZeros --centerReads --ignoreDuplicates --minMappingQuality 10"
    wrapper:
        "v1.31.1/bio/deeptools/plotcoverage"
    