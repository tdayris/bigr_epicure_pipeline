rule chipseeker_plot_annobar_single_sample:
    input:
        ranges="chipseeker/annotation/{sample}.{peaktype}.RDS",
    output:
        png="data_output/Peak_Calling/{peaktype}/Feature_Distribution/{sample}.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annobar/{sample}.{peaktype}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_annobar.R"


rule chipseeker_plot_annobar_differential_binding:
    input:
        ranges="chipseeker/annotation/{comparison_name}.RDS",
    output:
        png="data_output/Differential_Binding/{comparison_name}/Feature_Distribution.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annobar/{comparison_name}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_annobar.R"
