rule samtools_cram:
    input:
        "sambamba/markdup/{sample}.bam",
        bam_idx="sambamba/markdup/{sample}.bam.bai",
        ref=genome_fasta_path,
        ref_fai=genome_fai_path,
    output:
        protected("data_output/CRAM/{sample}.cram"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "samtools/view/{sample}.cram.log",
    params:
        extra="",
        region="",
    wrapper:
        "v1.31.1/bio/samtools/view"
