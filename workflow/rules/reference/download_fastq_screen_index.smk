rule download_fastq_screen_index:
    output:
        directory("reference/fastq_screen/index/{fq_genome}"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    cache: "omit-software"
    log:
        "logs/fastq_screen/download/{fq_genome}.log",
    params:
        extra=(
            "--no-check-certificate "
            "--recursive --no-parent "
            "--reject 'index.html*'"
        ),
        url="http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes/{fq_genome}",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params} --directory-prefix {output} "
        "{params.url} > {log} 2>&1"


rule download_fastq_screen_index_bisulfite:
    output:
        directory("reference/fastq_screen/index_bs/{fq_genome}"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    cache: "omit-software"
    log:
        "logs/fastq_screen/download/{fq_genome}.log",
    params:
        extra=(
            "--no-check-certificate "
            "--recursive --no-parent "
            "--reject 'index.html*'"
        ),
        url="http://ftp1.babraham.ac.uk/ftpusr46/FastQ_Screen_Genomes_Bisulfite/{fq_genome}",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params} --directory-prefix {output} "
        "{params.url} > {log} 2>&1"
