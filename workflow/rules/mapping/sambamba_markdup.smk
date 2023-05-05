rule sambamba_markdup:
    input:
        "samtools/view/{sample}.bam",
    output:
        temp("sambamba/markdup/{sample}.bam"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra="--remove-duplicates --overflow-list-size 600000",
    log:
        "logs/sambamba/markdup/{sample}.log",
    wrapper:
        "v1.29.0/bio/sambamba/markdup"
