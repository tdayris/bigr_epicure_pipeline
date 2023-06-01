# WARNING: release >= 75
rule blacklist_grch38:
    output:
        "reference/blacklist/homo_sapiens.GRCh38.{release}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh38.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg38-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: 67 < release <= 102
rule blacklist_mm10:
    output:
        "reference/blacklist/mus_musculus.GRCm38.102.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/mus_musculus.GRCm38.102.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm10-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: release <= 75
rule blacklist_grch37:
    output:
        "reference/blacklist/homo_sapiens.GRCh37.{release}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh37.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg19-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


# WARNING: release should be <= 67
rule blacklist_mm9:
    output:
        "reference/blacklist/mus_musculus.NCBIM37.{release}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/mus_musculus.NCBIM37.{release}.log",
    cache: "omit-software"
    params:
        address="https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
        extra="--verbose",
    conda:
        "../../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} -O {output} > {log} 2>&1"


rule bedtools_merge_blacklist:
    input:
        "reference/blacklist/{species}.{build}.{release}.bed.gz",
    output:
        "reference/blacklist/{species}.{build}.{release}.merged.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir=tmp,
    log:
        "logs/bedtools/merge/blacklist/{species}.{build}.{release}.log",
    params:
        extra="-d 5",
    wrapper:
        "v1.31.1/bio/bedtools/merge"
