rule bedtools_intersect_macs2:
    input:
        unpack(get_bedtools_intersect_macs2_input),
    output:
        temp("bedtools/intersect/{sample}.{peaktype}.bam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/bedtools/intersect/{sample}.{peaktype}.log",
    params:
        extra=" -wa ",
    wrapper:
        "v1.32.1/bio/bedtools/intersect"
