rule chipseeker_plot_distance_to_tss_single_sample:
    input:
        ranges="chipseeker/annotation/{sample}.{peaktype}.RDS",
    output:
        png="data_output/PeakCalling/{peaktype}/Distance_to_TSS/{sample}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annotate/{sample}.{peaktype}.log"
    params:
        extra=""
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot.R"


rule chipseeker_plot_distance_to_tss_differential_binding:
    input:
        ranges="chipseeker/annotation/{comparison_name}.RDS",
    output:
        png="data_output/DifferentialBinding/{comparison_name}/Distance_to_TSS.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annotate/{comparison_name}.log"
    params:
        extra=""
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot.R"