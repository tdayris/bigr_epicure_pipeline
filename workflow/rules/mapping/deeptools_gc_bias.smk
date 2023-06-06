rule deeptools_estimate_gc_bias:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bai="sambamba/markdup/{sample}.bam.bai",
        blacklist=blacklist_path,
        genome=genome_twobit_path,
    output:
        freq=temp("deeptools/gc_bias/{sample}.freq.txt"),
        plot="data_output/QC/GC_Bias/{sample}.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/compute_gc_bias/{sample}.log",
    params:
        extra="--plotFileFormat png ",
        genome_size=effective_genome_size,
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "computeGCBias "
        "--bamfile {input.bam} "
        "--genome {input.genome} "
        "--GCbiasFrequenciesFile {output.freq} "
        "--effectiveGenomeSize {params.genome_size} "
        "--blackListFileName {input.blacklist} "
        "--numberOfProcessors {threads} "
        "--biasPlot {output.plot} "
        "{params.extra} "
        "> {log} 2>&1 "


rule correct_gc_bias:
    input:
        bam="sambamba/markdup/{sample}.bam",
        bai="sambamba/markdup/{sample}.bam.bai",
        genome=genome_twobit_path,
        freq="deeptools/gc_bias/{sample}.freq.txt",
    output:
        temp("deeptools/corrected/{sample}.bam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/deeptools/correct_gc_bias/{sample}.log",
    params:
        extra="",
        genome_size=effective_genome_size,
    conda:
        "../../envs/deeptools.yaml"
    shell:
        "correctGCBias "
        "--bamfile {input.bam} "
        "--genome {input.genome} "
        "--GCbiasFrequenciesFile {input.freq} "
        "--correctedFile {output} "
        "--effectiveGenomeSize {params.genome_size} "
        "--numberOfProcessors {threads} "
        "{params.extra} "
        "> {log} 2>&1 "
