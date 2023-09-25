rule chipseeker_genome_cov_single_sample:
    input:
        unpack(get_chipseeker_genome_cov_single_sample_input),
    output:
        png="data_output/Peak_Calling/{peaktype}/Genome_Coverage/{sample}.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genomecov/{sample}.{peaktype}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_covplot.R"


rule chipseeker_plot_genome_cov_differential_binding:
    input:
        unpack(get_chipseeker_annotate_peak_from_ranges_input),
    output:
        png="data_output/Differential_Binding/{model_name}/Genome_Coverage.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/chipseeker/genomecov/{model_name}.log",
    params:
        extra="",
    conda:
        "../../envs/chipseeker.yaml"
    script:
        "../../scripts/chipseeker/chipseeker_covplot.R"
