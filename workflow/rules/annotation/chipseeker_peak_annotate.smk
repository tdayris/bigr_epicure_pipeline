rule chipseeker_annotate_peak:
    input:
        ranges="csaw/results/{comparison_name}.RDS"
    output:
        rds="chipseeker/annotation/{comparison_name}.RDS",
        tsv="data_output/Differential_Binding/{comparison_name}.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annotate/{comparison_name}.log"
    params:
        organism="hg38" if config[]