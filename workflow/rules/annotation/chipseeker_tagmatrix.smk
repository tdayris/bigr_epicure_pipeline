rule chipseeker_tagmatrix_from_ranges:
    input:
        ranges="csaw/results/{model_name}.RDS",
    output:
        rds="chipseeker/tagmatrix/{model_name}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 40,
        tmpdir=tmp,
    log:
        "logs/chipseeker/tagmatrix/{model_name}.log",
    params:
        organism="hg38" if species == "homo_sapiens" else "mm10",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_make_tagmatrix.R"


rule chipseeker_tagmatrix_from_bed:
    input:
        bed="macs2/callpeak_{peaktype}/{sample}_peaks.{peaktype}Peak.bed",
    output:
        rds="chipseeker/tagmatrix/{sample}.{peaktype}.RDS",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 40,
        tmpdir=tmp,
    log:
        "logs/chipseeker/tagmatrix/{sample}.{peaktype}.log",
    params:
        organism="hg38" if species == "homo_sapiens" else "mm10",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_make_tagmatrix.R"
