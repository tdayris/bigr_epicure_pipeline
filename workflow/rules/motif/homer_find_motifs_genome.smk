rule add_chr_bed:
    input:
        unpack(get_peak_file_list),
    output:
        temp("homer/peaks/{sample}.{peaktype}.bed"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 1,
        tmpdir=tmp,
    log:
        "logs/homer/add_chr/{sample}.{peaktype}.log",
    params:
        script=workflow.source_path("../../scripts/misc/bed_add_chr_to_sequence.awk"),
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.script} {input} > {output} 2> {log}"


rule homer_find_motif_genome:
    input:
        unpack(get_homer_find_motif_input),
    output:
        motif="homer/motif/{peaktype}/{sample}/homerMotifs.motifs",
        all_motif="homer/motif/{peaktype}/{sample}/homerMotifs.all.motifs",
        params=temp("homer/motif/{peaktype}/{sample}/motifFindingParameters.txt"),
        known="homer/motif/{peaktype}/{sample}/knownResults.txt",
        seq="homer/motif/{peaktype}/{sample}/seq.autonorm.tsv",
        html="data_output/Motifs/{peaktype}/{sample}/homerResults.html",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 60 * 2,
        tmpdir=tmp,
    log:
        "logs/homer/find_motif_genome/{sample}.{peaktype}.log",
    params:
        genome="hg38" if build == "GRCh38" else "mm10",
        size=200,
        output_directory="data_output/Motifs/{peaktype}/{sample}",
    conda:
        "../../envs/homer.yaml"
    shell:
        "findMotifsGenome.pl "
        "{input.peak} "
        "{params.genome} "
        "{params.output_directory} "
        "-size {params.size} "
        "> {log} 2>&1 "
