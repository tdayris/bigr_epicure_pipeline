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
        "v1.29.0/bio/samtools/view"


rule samtools_cram_shifted:
    input:
        "deeptools/alignment_sieve/{sample}.bam",
        bam_idx="deeptools/alignment_sieve/{sample}.bam.bai",
        ref=genome_fasta_path,
        ref_fai=genome_fai_path,
    output:
        protected("data_output/CRAM-shifted/{sample}.cram"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "samtools/view/{sample}.shifted.cram.log",
    params:
        extra="",
        region="",
    wrapper:
        "v1.29.0/bio/samtools/view"
