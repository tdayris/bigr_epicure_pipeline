rule chipseeker_genebody_cov_single_sample:
    input:
        ranges="chipseeker/tagmatrix/{sample}.{peaktype}.RDS",
    output:
        png="data_output/Peak_Calling/{peaktype}/Gene_Body_Coverage/{sample}.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genebody_cov/{sample}.{peaktype}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_gene_body.R"


rule chipseeker_genebody_cov_differential_binding:
    input:
        ranges="chipseeker/annotation/{model_name}.RDS",
    output:
        png="data_output/Differential_Binding/{model_name}/Gene_Body_Coverage.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genebody_cov/{model_name}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_plot_gene_body.R"
