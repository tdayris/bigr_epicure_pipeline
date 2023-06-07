rule homer_annotatepeaks:
    input:
        peaks="macs2/callpeak_{peaktype}/{sample}_peaks.{peaktype}Peak.bed",
        genome=genome_fasta_path,
        motif_files="homer/motif/{peaktype}/{sample}/homerMotifs.motifs",
        gtf=genome_annotation_path,
        wig="data_output/Coverage/{sample}.bw",
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
    log:
        "logs/homer/annotate/{sample}.{peaktype}.log",
    wrapper:
        "v1.31.1/bio/homer/annotatePeaks"
