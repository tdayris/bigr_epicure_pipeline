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
        directory("data_output/QC/Trimming.QC_data"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    params:
        extra=" --module fastp --module fastq_screen ",
        use_input_files_only=True,
    log:
        "logs/multiqc/trimming.log",
    wrapper:
        f"{snakemake_wrappers_version}/bio/multiqc"
