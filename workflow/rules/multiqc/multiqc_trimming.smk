rule multiqc_trimming:
    input:
        unpack(get_multiqc_trimming_input),
    output:
        report(
            "data_output/QC/Trimming.QC.html",
            caption="../../report/multiqc/trimming.rst",
            category="Quality Control",
            labels={"type": "html", "category": "QC"},
        ),
        directory("data_output/QC/Trimming.QC.data"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra=" --module fastp --module fastq_screen ",
    log:
        "logs/multiqc/trimming.log",
    wrapper:
        "v1.31.1/bio/multiqc"
