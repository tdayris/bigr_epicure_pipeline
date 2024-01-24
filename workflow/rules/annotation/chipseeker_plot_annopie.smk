rule chipseeker_plot_annopie_single_sample:
    input:
        ranges="chipseeker/annotation/{sample}.{peaktype}.RDS",
    output:
        png=report(
            "data_output/Peak_Calling/{peaktype}/Feature_Distribution_Pie/{sample}.png",
            caption="../../report/annotation/genome_annot_single_sample.rst",
            category="Annotation",
            labels={"type": "png", "category": "Annotation"},
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annopie/{sample}.{peaktype}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_annopie.R"


rule chipseeker_plot_annopie_differential_binding:
    input:
        ranges="chipseeker/annotation/{model_name}.RDS",
    output:
        png="data_output/Differential_Binding/{model_name}/Feature_Distribution_Pie.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annopie/{model_name}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_annopie.R"


rule chipseeker_plot_annopie_single_sample_list:
    input:
        ranges=expand(
            "chipseeker/annotation/{sample}.{peaktype}.RDS",
            sample=get_tested_sample_list(),
            allow_missing=True,
        ),
    output:
        png="data_output/Peak_Calling/{peaktype}/GenomicAnnotations_Pie.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/annopie/all.{peaktype}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_annopie_list.R"
