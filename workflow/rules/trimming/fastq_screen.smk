rule fastq_screen_paired:
    input:
        unpack(get_fastq_screen_input),
    output:
        png=temp("data_output/{sample}.{stream}.png"),
        txt=temp("fastq_screen/{sample}.{stream}.txt"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 12,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fastq_screen/{sample}.{stream}.log",
    params:
        fastq_screen_config=lambda wildcards, input: str(input[1]),
        subset=100000,
        aligner="bowtie2",
    wrapper:
        "1.29.0/bio/fastq_screen"


rule fastq_screen_single:
    input:
        unpack(get_fastq_screen_input),
    output:
        png=temp("data_output/{sample}.png"),
        txt=temp("fastq_screen/{sample}.txt"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 12,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/fastq_screen/{sample}.log",
    params:
        fastq_screen_config=lambda wildcards, input: str(input[1]),
        subset=100000,
        aligner="bowtie2",
    wrapper:
        "1.29.0/bio/fastq_screen"
