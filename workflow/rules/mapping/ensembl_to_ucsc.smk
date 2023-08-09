rule sambamba_view_unzip_bam_header:
    input:
        "sambamba/markdup/{sample}.bam"
    output:
        temp("sambamba/markdup/{sample}.header")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/unzip/{sample}.header.log",
    group:
        "reheader_as_ucsc"
    params:
        extra="-H",
    wrapper:
        "v2.0.0/bio/sambamba/view"


rule sed_edit_header_sequences:
    input:
        "sambamba/markdup/{sample}.header"
    output:
        temp("sambamba/markdup/{sample}.sq.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 55,
        tmpdir=tmp,
    group:
        "reheader_as_ucsc"
    log:
        "logs/sed/header/{sample}.log"
    params:
        "'s/SN:/SN:chr/g'"
    conda:
        "../../envs/bash.yaml"
    shell:
        "sed {params} {input} > {output} 2> {log}"


rule sambamba_view_unzip_bam:
    input:
        "sambamba/markdup/{sample}.bam"
    output:
        temp("sambamba/markdup/{sample}.content.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/unzip/{sample}.content.log",
    group:
        "ensembl_as_ucsc"
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/view"


rule perl_edit_read_sequences:
    input:
        "sambamba/markdup/{sample}.content.sam"
    output:
        temp("sambamba/markdup/{sample}.ucsc.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/sequences/{sample}.log"
    params:
        awk=workflow.source_path("../../scripts/formatter/ensembl_to_ucsc.awk"),
    group:
        "ensembl_as_ucsc"
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.awk} {input} > {output} 2> {log}"


rule cat_sam_file:
    input:
        header="sambamba/markdup/{sample}.sq.sam",
        reads="sambamba/markdup/{sample}.ucsc.sam",
    output:
        temp("sambamba/ucsc/{sample}.sam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime= lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/cat/ucsc/{sample}.log"
    params:
        extra=""
    conda:
        "../../envs/bash.yaml"
    group:
        "compress_ucsc"
    shell:
        "cat {params.extra} {input.header} {input.reads} > {output} 2> {log}"

rule sambamba_zip_sam:
    input:
        "sambamba/ucsc/{sample}.sam",
    output:
        temp("sambamba/ucsc/{sample}.unsorted.bam"),
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 35,    
        tmpdir=tmp,
    group:
        "compress_ucsc"
    log:
        "logs/sambamba/ucsc/{sample}.zip.log"
    params:
        extra="--format bam --with-header"
    wrapper:
        "v2.0.0/bio/sambamba/view"

rule sambamba_sort_ucsc:
    input:
        "sambamba/ucsc/{sample}.unsorted.bam"
    output:
        "sambamba/ucsc/{sample}.bam"
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{sample}.ucsc.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/sort"
rule sambamba_index_ucsc:
    input:
        "sambamba/ucsc/{sample}.bam"
    output:
        "sambamba/ucsc/{sample}.bam.bai"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/sambamba/ucsc/{sample}.index.log"
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/index"


rule sambamba_view_unzip_bam_headerfactor:
    input:
        "sambamba/markdup/{factor_level}.bam"
    output:
        temp("sambamba/markdup/{factor_level}.header")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/unzip/{factor_level}.header.log",
    group:
        "reheader_as_ucsc"
    params:
        extra="-H",
    wrapper:
        "v2.0.0/bio/sambamba/view"


rule sed_edit_header_sequences_factor_level:
    input:
        "sambamba/markdup/{factor_level}.header"
    output:
        temp("sambamba/markdup/{factor_level}.sq.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 55,
        tmpdir=tmp,
    group:
        "reheader_as_ucsc"
    log:
        "logs/sed/header/{factor_level}.log"
    params:
        "'s/SN:/SN:chr/g'"
    conda:
        "../../envs/bash.yaml"
    shell:
        "sed {params} {input} > {output} 2> {log}"


rule sambamba_view_unzip_bam_factor:
    input:
        "sambamba/markdup/{factor_level}.bam"
    output:
        temp("sambamba/markdup/{factor_level}.content.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 20,
        tmpdir=tmp,
    log:
        "logs/sambamba/view/unzip/{factor_level}.content.log",
    group:
        "ensembl_as_ucsc"
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/view"



rule perl_edit_read_sequences_per_factor:
    input:
        "sambamba/markdup/{factor_level}.content.sam"
    output:
        temp("sambamba/markdup/{factor_level}.ucsc.sam")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/sequences/{factor_level}.log"
    params:
        awk=workflow.source_path("../../scripts/formatter/ensembl_to_ucsc.awk"),
    group:
        "ensembl_as_ucsc"
    conda:
        "../../envs/bash.yaml"
    shell:
        "awk --file {params.awk} {input} > {output} 2> {log}"

rule cat_sam_file_factor:
    input:
        header="sambamba/markdup/{factor_level}.sq.sam",
        reads="sambamba/markdup/{factor_level}.ucsc.sam",
    output:
        temp("sambamba/ucsc/{factor_level}.sam"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime= lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/cat/ucsc/{factor_level}.log"
    params:
        extra=""
    conda:
        "../../envs/bash.yaml"
    group:
        "compress_ucsc"
    shell:
        "cat {params.extra} {input.header} {input.reads} > {output} 2> {log}"

rule sambamba_zip_sam_factor:
    input:
        "sambamba/ucsc/{factor_level}.sam",
    output:
        temp("sambamba/ucsc/{factor_level}.unsorted.bam"),
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 35,    
        tmpdir=tmp,
    group:
        "compress_ucsc"
    log:
        "logs/sambamba/ucsc/{factor_level}.zip.log"
    params:
        extra="--format bam --with-header"
    wrapper:
        "v2.0.0/bio/sambamba/view"

rule sambamba_sort_ucsc:
    input:
        "sambamba/ucsc/{factor_level}.unsorted.bam",
    output:
        "sambamba/ucsc/{factor_level}.bam"
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20 * 1024,
        runtime=lambda wildcards, attempt: attempt * 60 * 5,
        tmpdir=tmp,
    log:
        "logs/sambamba/sort/{factor_level}.ucsc.log",
    params:
        extra=lambda wildcards, resources: f"--memory-limit {resources.mem_mb - 1024}MiB",
    wrapper:
        f"{snakemake_wrappers_version}/bio/sambamba/sort"

rule sambamba_index_ucsc:
    input:
        "sambamba/ucsc/{factor_level}.bam"
    output:
        "sambamba/ucsc/{factor_level}.bam.bai"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/sambamba/ucsc/{factor_level}.index.log"
    params:
        extra="",
    wrapper:
        "v2.0.0/bio/sambamba/index"
