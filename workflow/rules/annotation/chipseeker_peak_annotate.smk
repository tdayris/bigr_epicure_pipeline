rule chipseeker_annotate_peak_from_ranges:
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
        organism="hg38" if species == "homo_sapiens" else "mm10"
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_annotate.R"


rule chipseeker_annotate_peak_from_bed:
    input:
        bed="macs2/callpeak_{peaktype}/{sample}_peaks.{peaktype}Peak.bed"
    output:
        rds="chipseeker/annotation/{sample}.{peaktype}.RDS",
        tsv="data_output/Differential_Binding/{sample}.{peaktype}.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annotate/{sample}.{peaktype}.log"
    params:
        organism="hg38" if species == "homo_sapiens" else "mm10"
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_annotate.R"