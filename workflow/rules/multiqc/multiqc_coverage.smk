rule multiqc_coverage:
    input:
        unpack(get_multiqc_coverage_input),
    output:
        report(
            "data_output/QC/Coverage.QC.html",
            caption="../../report/multiqc/coverage.rst",
            category="Quality Control",
            labels={"type": "html", "category": "QC"},
        ),
        directory("data_output/QC/Coverage.QC_data"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra=" --module fastp --module fastq_screen --module samtools --module picard --module deeptools",
        use_input_files_only=True,
    log:
        "logs/multiqc/coverage.log",
    wrapper:
        "v1.31.1/bio/multiqc"
