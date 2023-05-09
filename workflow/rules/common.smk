import logging
import os
import pandas

from snakemake.remote import FTP
from typing import Any, Dict


########################
### Load design file ###
########################


design = pandas.read_csv(
    config.get("design", "../config/design.tsv"),
    sep="\t",
    header=0,
    index_col=0,
)


########################
### Global variables ###
########################

# FTP handler for possible downloads
FTP = FTP.RemoteProvider()


# This let all temporary file end-up in the same directory
tmp = os.environ.get("TMP", "temp")


# Mouse has less chromosomes, but this will work anyway
# Does not work for organisms with more than 22 autosomal chromosomes
# Does not work for organisms with no X/Y system for sexual chromosomes
canonical_chromosomes += list(map(str, range(1, 23))) + ["X", "Y"]
canonical_chromosomes += [f"chr{chrom}" for chrom in canonical_chromosomes]


# Protocol. It can be either : `chip-seq`, `atac-seq`, `cut&run`, or `cut&tag`
protocol = config.get("protocol", "chip-seq")

################################
### Paths to reference files ###
################################

# Main genome informations
species = config.get("reference", {}).get("species", "homo_sapiens")
build = config.get("reference", {}).get("build", "GRCh38")
release = config.get("reference", {}).get("release", "109")

# Genome sequence (FASTA, not gzipped)
genome_fasta_path = config.get("reference", {}).get("genome_fasta")
if not genome_fasta_path:
    logging.info(
        "Missing genome FASTA path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_fasta_path = f"reference/{species}.{build}.{release}.fasta"


# Genome sequences indexes
genome_fai_path = genome_fasta_path + ".fai"
genome_dict_path = ".".join(genome_fasta_path.split(".")[:-1]) + ".dict"


# Genome annotation (GFF/GTF, not gzipped)
genome_annotation_path = config.get("reference", {}).get("genome_gtf")
if not genome_annotation_path:
    logging.info(
        "Missing GTF path in the file `config.yaml`. "
        "A new one will be downloaded from Ensembl FTP."
    )
    genome_annotation_path = f"reference/{species}.{build}.{release}.gtf"


# Bowtie2 genome index
bowtie2_index_path = config.get("reference", {}).get("bowtie2_index")
if not bowtie2_index_path:
    logging.info(
        "Missing Bowtie2 index path in the file `config.yaml`. "
        "A new one will be created."
    )
    bowtie2_index = multiext(
        f"reference/bowtie2_index/{species}.{build}.{release}",
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )


# Genome blacklisted regions
blacklist_path = config.get("reference", {}).get("blacklist")
if not blacklist_path:
    blacklist_path = f"reference/blacklist/{species}.{build}.{release}.bed.gz"


# See: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
default_effective_genome_size = {
    "GRCz10": {
        "50": 1_195_445_591,
        "75": 1_251_132_686,
        "100": 1_280_189_044,
        "150": 1_312_207_169,
        "200": 1_321_355_241,
    },
    "WBcel235": {
        "50": 95_159_452,
        "75": 96_945_445,
        "100": 98_259_998,
        "150": 98_721_253,
        "200": 98_672_758,
    },
    "dm3": {
        "50": 130_428_560,
        "75": 135_004_462,
        "100": 139_647_232,
        "150": 144_307_808,
        "200": 148_524_010,
    },
    "dm6": {
        "50": 125_464_728,
        "75": 127_324_632,
        "100": 129_789_873,
        "150": 129_941_135,
        "200": 132_509_163,
    },
    "GRCh37": {
        "50": 2_685_511_504,
        "75": 2_736_124_973,
        "100": 2_776_919_808,
        "150": 2_827_437_033,
        "200": 2_855_464_000,
    },
    "GRCh38": {
        "50": 2_701_495_761,
        "75": 2_747_877_777,
        "100": 2_805_636_331,
        "150": 2_862_010_578,
        "200": 2_887_553_303,
    },
    "GRCm37": {
        "50": 2_304_947_926,
        "75": 2_404_646_224,
        "100": 2_462_481_010,
        "150": 2_489_384_235,
        "200": 2_513_019_276,
    },
    "GRCm38": {
        "50": 2_308_125_349,
        "75": 2_407_883_318,
        "100": 2_467_481_108,
        "150": 2_494_787_188,
        "200": 2_520_869_189,
    },
}
read_length = config.get("library", {}).get("read_length", 100)
effective_genome_size = default_effective_genome_size.get(build).get(
    read_length, 2_805_636_331
)

###################
### IO function ###
###################


def get_fastp_input(
    wildcards, design: pandas.DataFrame = design
) -> Dict[str, List[str]]:
    """
    Return the list of Fastp input files
    """
    fastp_input = {
        "sample": [design["Upstream_file"].loc[wildcards.sample]],
    }
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if down is not None and down != "":
            fastp_input["sample"].append(down)
    return fastp_input


def get_bowtie2_align_input(
    wildcards, design: pandas.DataFrame = design, bowtie2_index: str = bowtie2_index
) -> Dict[str, List[str]]:
    """
    Return the list of bowtie2 align input
    """
    bowtie2_align_input = {"idx": bowtie2_index_path}
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if (down is not None) and (down != ""):
            bowtie2_align_input["sample"] = expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                sample=[wildcards.sample],
            )
        else:
            bowtie2_align_input["sample"] = expand(
                "fastp/trimmed/se/{sample}.fastq", sample=[wildcards.sample]
            )

    return bowtie2_align_input


def get_fastp_output_html(design: pandas.DataFrame = design) -> List[str]:
    """
    Return the list of fastp reports
    """
    html_list = []
    for sample in design.index:
        down = design["Downstream_file"].loc[sample]
        if (down is not None) and (down != ""):
            html_list.append("data_output/qc/fastp/{sample}.se.html")
        else:
            html_list.append("data_output/qc/fastp/{sample}.pe.html")
    return html_list


def get_deeptools_bamcoverage_input(
    wildcards, protocol: str = protocol, blacklist_path: str = blacklist_path
) -> Dict[str, str]:
    """
    Return the expected list of DeepTools input file list
    """
    deeptools_bamcoverage_input = {"blacklist": blacklist_path}
    if protocol == "ataq-seq":
        deeptools_bamcoverage_input["bam"] = "deeptools/alignment_sieve/{sample}.bam"
        deeptools_bamcoverage_input[
            "bai"
        ] = "deeptools/alignment_sieve/{sample}.bam.bai"
    else:
        deeptools_bamcoverage_input["bam"] = "sambamba/markdup/{sample}.bam"
        deeptools_bamcoverage_input["bai"] = "sambamba/markdup/{sample}.bam.bai"

    return deeptools_bamcoverage_input


def get_fastq_screen_input(
    wildcards,
    config: Dict[str, Any] = config,
    design: pandas.DataFrame = design,
    default_fastq_screen_genomes: List[str] = default_fastq_screen_genomes,
    fastq_screen_index_path: str = fastq_screen_index_path,
) -> List[str]:
    """
    Return the list of expected input files for fastq_screen

    First file in the input file list MUST be the fastq file of interest.
    Second file in the input file list MUST be the configuration file.
    """
    fastq_screen_input = []

    # Let the fastq file be the first file in the input file list
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if (down is not None) and (down != ""):
            fastq_screen_input.append(
                f"fastp/trimmed/pe/{wildcards.sample}.{wildcards.stream}.fastq"
            )
        else:
            fastq_screen_input.append(f"fastp/trimmed/se/{wildcards.sample}.fastq")

    # Let the config file be the second input
    fq_conf = config.get("reference", {}).get(
        "fastq_screen_config", "../../../config/fastq_screen.conf"
    )
    fastq_screen_input.append("../../../config/fastq_screen.conf")

    # Optional other input files
    if config.get("reference", {}).get("download_fastq_screen_indexes", False):
        for fq_genome in fastq_screen_genomes:
            fastq_screen_input.append(
                f"reference/fastq_screen/index/{fq_genome}",
            )

    return fastq_screen_input


def get_macs2_callpeak_input(wildcards, protocol: str = "chip-seq") -> Dict[str, str]:
    """
    Return expected list of input files for Macs2 callpeak
    """
    macs2_callpeak_input = {}
    if protocol == "atac-seq":
        macs2_callpeak_input["teatment"] = "deeptools/alignment_sieve/{sample}.bam"
    else:
        macs2_callpeak_input["teatment"] = "sambamba/markdup/{sample}.bam"

    if "Input" in design.keys():
        down = design["Input"].loc[wildcards.sample]
        if (down is not None) and (down != ""):
            if protocol == "atac-seq":
                macs2_callpeak_input["control"] = "deeptools/alignment_sieve/{down}.bam"
            else:
                macs2_callpeak_input["control"] = "sambamba/markdup/{down}.bam"

    return macs2_callpeak_input


def get_macs2_params(
    wildcards, effective_genome_size: int = effective_genome_size
) -> str:
    """
    Return expected parameters for Macs2 callpeak
    """
    extra = f" --gsize {effective_genome_size} "
    if "Downstream_file" in design.keys():
        down = design["Downstream_file"].loc[wildcards.sample]
        if (down is not None) and (down != ""):
            extra += " --format BAMPE "
    else:
        extra += " --format BAM "

        if "Fragment_size" in design.keys():
            fs = design["Fragment_size"].loc[wildcards.sample]
            if (down is not None) and (down != ""):
                extra += f" --nomodel --extsize {fs} "
            else:
                raise ValueError(
                    "Single-ended reads should have a "
                    "`Fragment_size` associated in the design file."
                )

    return extra


def targets(
    config: Dict[str, Any] = config,
    design: pandas.DataFrame = design,
    protocol: str = protocol,
    genome_fasta_path: str = genome_fasta_path,
    genome_annotation_path: str = genome_annotation_path,
    bowtie2_index_path: str = bowtie2_index_path,
):
    """
    Return the list of expected output files, depending on the
    choices made by user in configuration file at: `<root>/config/config.yaml`
    """
    expected_targets = {}
    steps = config.get("steps", {})

    if steps.get("install", False):
        expected_targets["gtf"] = genome_annotation_path
        expected_targets["genome_fasta"] = genome_fasta_path
        expected_targets["bowtie2_index"] = bowtie2_index_path
        expected_targets["fastq_screen_indexes"] = (
            expand(
                "reference/fastq_screen/index/{fq_genome}",
                fq_genome=fastq_screen_genomes,
            ),
        )

    if steps.get("trimming", False):
        expected_targets["fastp"] = get_fastp_output_html(design=design)

    if steps.get("mapping", False):
        if protocol == "atac-seq":
            expected_targets["mapping"] = "data_output/CRAM-shifted/{sample}.cram"
        else:
            expected_targets["mapping"] = "data_output/CRAM/{sample}.cram"

    if steps.get("coverage", False):
        raise NotImplementedError("Coverage analysis not yet implements")

    if steps.get("calling", False):
        raise NotImplementedError("Peak calling not yet implements")

    if steps.get("diff_cov", False):
        raise NotImplementedError("Differential coverage analysis not yet implements")

    if steps.get("motives", False):
        raise NotImplementedError("Mitives analysis not yet implements")
