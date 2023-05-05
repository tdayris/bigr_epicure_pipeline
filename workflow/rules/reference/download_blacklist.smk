# WARNING: build >= 75
rule blacklist_grch38:
    input:
        FTP.remote(
            "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg38-blacklist.bed.gz",
            static=True,
            keep_local=True,
        ),
    output:
        "reference/blacklist/homo_sapiens.GRCh38.{build}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh38.{build}.log",
    cache: "omit-software"
    params:
        "--verbose",
    shell:
        "bash.yaml"
    shell:
        "mv {params} {input} {output} > {log} 2>&1"


# WARNING: 67 < build <= 102
rule blacklist_mm10:
    input:
        FTP.remote(
            "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm10-blacklist.bed.gz",
            static=True,
            keep_local=True,
        ),
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
        "--verbose",
    shell:
        "bash.yaml"
    shell:
        "mv {params} {input} {output} > {log} 2>&1"


# WARNING: build <= 75
rule blacklist_grch37:
    input:
        FTP.remote(
            "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/hg19-blacklist.bed.gz",
            static=True,
            keep_local=True,
        ),
    output:
        "reference/blacklist/homo_sapiens.GRCh37.{build}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/homo_sapiens.GRCh37.{build}.log",
    cache: "omit-software"
    params:
        "--verbose",
    shell:
        "bash.yaml"
    shell:
        "mv {params} {input} {output} > {log} 2>&1"


# WARNING: build should be <= 67
rule blacklist_mm9:
    input:
        FTP.remote(
            "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz",
            static=True,
            keep_local=True,
        ),
    output:
        "reference/blacklist/mus_musculus.NCBIM37.{build}.bed.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 60 * 24,
        tmpdir=tmp,
    log:
        "logs/ftp/blacklist/mus_musculus.NCBIM37.{build}.log",
    cache: "omit-software"
    params:
        "--verbose",
    shell:
        "bash.yaml"
    shell:
        "mv {params} {input} {output} > {log} 2>&1"
