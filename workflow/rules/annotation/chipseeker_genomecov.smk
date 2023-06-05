rule chipseeker_genome_cov_single_sample:
    input:
        ranges="chipseeker/annotation/{sample}.{peaktype}.RDS",
    output:
        png="data_output/Peak_Calling/{peaktype}/Genome_Coverage/{sample}.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genomecov/{sample}.{peaktype}.log"
    params:
        extra=""
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_covplot.R"


rule chipseeker_plot_genome_cov_differential_binding:
    input:
        ranges="chipseeker/annotation/{comparison_name}.RDS",
    output:
        png="data_output/Differential_Binding/{comparison_name}/Genome_Coverage.png"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genomecov/{comparison_name}.log"
    params:
        extra=""
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_covplot.R"