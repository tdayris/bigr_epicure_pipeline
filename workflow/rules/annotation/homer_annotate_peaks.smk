rule homer_annotatepeaks:
    input:
        unpack(get_homer_annotatepeaks_input),
    output:
        annotations="data_output/Motifs/{peaktype}/{sample}_homer_annot.txt",
        matrix=multiext(
            "homet/{peaktype}/{sample}",
            ".count.matrix.txt",
            ".ratio.matrix.txt",
            ".logPvalue.matrix.txt",
            ".stats.txt",
        ),
        mfasta="data_output/Motifs/{peaktype}/{sample}_motif.fasta",
        mbed="data_output/Motifs/{peaktype}/{sample}_motif.bed",
        mlogic="data_output/Motifs/{peaktype}/{sample}_motif.logic",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    params:
        mode="tss hg38" if build == "GRCh38" else "tss mm10",
        extra="-CpG",
        genome="hg38" if build == "GRCh38" else "mm10",
    log:
        "logs/homer/annotate/{sample}.{peaktype}.log",
    conda:
        "../../envs/homer.yaml"
    shell:
        "annotatePeaks.pl "
        "{input.peak} "
        "{params.genome} "
        "-cpu {threads} "
        "-wig {input.wig} "
        "-m {input.motif_files} "
        "-mfasta {output.mfasta} "
        "-mbed {output.mbed} "
        "-mlogic {output.mlogic} "
        "> {log} 2>&1 "
